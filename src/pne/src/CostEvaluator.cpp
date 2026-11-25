// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/CostEvaluator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_set>
#include <vector>

#include "odb/dbTransform.h"
#include "utl/Logger.h"
#include "sta/Network.hh"

#include "pne/SoftMacro.h"

namespace pne {

namespace {

// Die-frame position of a cached macro pin, derived from the current
// B*-tree node coordinates and orientation.  The pin center is stored in
// master frame (R0); the same transform math as odb::dbInst placement is
// applied so results match dbITerm::getAvgXY() once the placement is
// committed to the DB.
odb::Point pinPositionFromNode(const BStarNode* node,
                               const odb::Point& master_pin,
                               int core_x,
                               int core_y)
{
  int mw = 0;
  int mh = 0;
  if (node->isSoftMacro()) {
    mw = node->getSoftMacro()->width;
    mh = node->getSoftMacro()->height;
  } else {
    odb::dbMaster* master = node->getInst()->getMaster();
    mw = static_cast<int>(master->getWidth());
    mh = static_cast<int>(master->getHeight());
  }

  odb::Point pin = master_pin;
  odb::Rect box(0, 0, mw, mh);
  const odb::dbTransform xform(node->getOrientation());
  xform.apply(pin);
  xform.apply(box);

  // Node (x, y) is the halo box lower-left; the macro bbox lower-left is
  // offset by the left/bottom halo.  Both are in tree-local coordinates,
  // so the placement core origin translates into the die frame.
  const Halo& halo = node->getHalo();
  const int x = core_x + node->getX() + halo.left + (pin.x() - box.xMin());
  const int y = core_y + node->getY() + halo.bottom + (pin.y() - box.yMin());
  return {x, y};
}

double overlapArea(const odb::Rect& a, const odb::Rect& b)
{
  const int64_t dx = static_cast<int64_t>(std::min(a.xMax(), b.xMax()))
                     - std::max(a.xMin(), b.xMin());
  const int64_t dy = static_cast<int64_t>(std::min(a.yMax(), b.yMax()))
                     - std::max(a.yMin(), b.yMin());
  return (dx > 0 && dy > 0) ? static_cast<double>(dx) * static_cast<double>(dy)
                            : 0.0;
}

// Bounding box of a point set; its HPWL is 0 while empty.
struct PointBox
{
  int64_t min_x = std::numeric_limits<int64_t>::max();
  int64_t max_x = std::numeric_limits<int64_t>::min();
  int64_t min_y = std::numeric_limits<int64_t>::max();
  int64_t max_y = std::numeric_limits<int64_t>::min();

  void add(int64_t x, int64_t y)
  {
    min_x = std::min(min_x, x);
    max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
  }

  int64_t hpwl() const
  {
    return min_x > max_x ? 0 : (max_x - min_x) + (max_y - min_y);
  }
};

}  // namespace

CostEvaluator::CostEvaluator(odb::dbDatabase* db,
                             sta::dbNetwork* network,
                             utl::Logger* logger)
    : db_(db), network_(network), logger_(logger)
{
}

void CostEvaluator::classifyNets(const std::vector<odb::dbInst*>& macros,
                                 BStarTree* tree)
{
  // Build macro lookup map
  macro_map_.clear();
  for (odb::dbInst* macro : macros) {
    macro_map_[macro] = true;
  }

  collectPlacementBlockages();

  // Map DB instances to tree nodes: hard macros map to their own node,
  // std cells inside a soft macro map to the soft macro's node (their pin
  // is approximated by the soft macro center).
  std::unordered_map<odb::dbInst*, int> inst_to_node;
  if (tree != nullptr) {
    for (const auto& node : tree->getNodes()) {
      if (node->isHardMacro()) {
        inst_to_node[node->getInst()] = node->getId();
      } else if (node->isSoftMacro()) {
        for (odb::dbInst* member : node->getSoftMacro()->instances) {
          inst_to_node[member] = node->getId();
        }
      }
    }
  }

  // Get block and iterate through all nets
  odb::dbBlock* block = db_->getChip()->getBlock();
  nets_.clear();
  net_to_idx_.clear();

  num_internal_nets_ = 0;
  num_io_nets_ = 0;
  num_hard_nets_ = 0;

  for (odb::dbNet* net : block->getNets()) {
    // Skip special nets (power, ground, clock)
    if (net->getSigType() == odb::dbSigType::POWER ||
        net->getSigType() == odb::dbSigType::GROUND ||
        net->getSigType() == odb::dbSigType::CLOCK) {
      continue;
    }

    NetInfo info;
    info.net = net;
    info.hpwl = 0;
    info.min_x = std::numeric_limits<int>::max();
    info.max_x = std::numeric_limits<int>::min();
    info.min_y = std::numeric_limits<int>::max();
    info.max_y = std::numeric_limits<int>::min();

    // Cache the net geometry: movable pins by tree node, everything else
    // folded into a fixed die-frame bounding box.  Unplaced instances
    // (std cells before global placement) carry no position information
    // and are skipped instead of polluting the bounding box with (0, 0).
    auto add_fixed_point = [&info](int x, int y) {
      if (!info.has_fixed) {
        info.has_fixed = true;
        info.fixed_min_x = info.fixed_max_x = x;
        info.fixed_min_y = info.fixed_max_y = y;
      } else {
        info.fixed_min_x = std::min(info.fixed_min_x, x);
        info.fixed_max_x = std::max(info.fixed_max_x, x);
        info.fixed_min_y = std::min(info.fixed_min_y, y);
        info.fixed_max_y = std::max(info.fixed_max_y, y);
      }
    };

    std::unordered_set<int> soft_nodes_seen;
    bool has_hard = false;
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();

      auto it = inst_to_node.find(inst);
      if (it != inst_to_node.end()) {
        const int node_id = it->second;
        const BStarNode* node = tree->getNodes()[node_id].get();
        if (node->isSoftMacro()) {
          // One point per soft macro per net is enough: all member pins
          // collapse to the region center.
          if (soft_nodes_seen.insert(node_id).second) {
            info.macro_pins.emplace_back(
                node_id,
                odb::Point(node->getSoftMacro()->width / 2,
                           node->getSoftMacro()->height / 2));
          }
        } else {
          odb::dbMTerm* mterm = iterm->getMTerm();
          odb::Rect pin_bbox = mterm->getBBox();
          info.macro_pins.emplace_back(
              node_id,
              odb::Point((pin_bbox.xMin() + pin_bbox.xMax()) / 2,
                         (pin_bbox.yMin() + pin_bbox.yMax()) / 2));
          has_hard = true;
        }
        continue;
      }

      if (!inst->isPlaced() && !inst->isFixed()) {
        continue;
      }

      int x, y;
      if (iterm->getAvgXY(&x, &y)) {
        add_fixed_point(x, y);
      } else {
        int inst_x, inst_y;
        inst->getLocation(inst_x, inst_y);
        add_fixed_point(inst_x, inst_y);
      }
    }

    // Nets that touch no movable macro are invariant under macro moves.
    if (info.macro_pins.empty()) {
      continue;
    }

    // Block terminals (chip IO pins) are fixed points at their current
    // (ppl-assigned) location.
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      for (odb::dbBPin* bpin : bterm->getBPins()) {
        for (odb::dbBox* box : bpin->getBoxes()) {
          odb::Rect rect = box->getBox();
          add_fixed_point((rect.xMin() + rect.xMax()) / 2,
                          (rect.yMin() + rect.yMax()) / 2);
        }
      }
    }

    info.type = classifyNet(net);
    if (has_hard && info.macro_pins.size() + (info.has_fixed ? 1 : 0) >= 2) {
      num_hard_nets_++;
    }

    net_to_idx_[net] = nets_.size();
    nets_.push_back(std::move(info));

    switch (nets_.back().type) {
      case NetType::INTERNAL:
        num_internal_nets_++;
        break;
      case NetType::IO:
        num_io_nets_++;
        break;
    }
  }

  logger_->info(utl::PNE, 10,
                "Net classification (macro-connected) - Internal: {}, IO: {}",
                num_internal_nets_, num_io_nets_);
}

NetType CostEvaluator::classifyNet(odb::dbNet* net)
{
  bool has_macro = false;
  bool has_io = false;
  bool has_std_cell = false;
  
  // Check instance terminals
  for (odb::dbITerm* iterm : net->getITerms()) {
    odb::dbInst* inst = iterm->getInst();
    const odb::dbMasterType master_type = inst->getMaster()->getType();
    
    if (macro_map_.find(inst) != macro_map_.end()) {
      has_macro = true;
    } else if (master_type.isPad()) {
      has_io = true;
    } else if (inst->isBlock()) {
      has_macro = true;  // Also consider blocks as macros
    } else {
      has_std_cell = true;
    }
  }
  
  // Check block terminals (IO pins)
  for (odb::dbBTerm* bterm : net->getBTerms()) {
    has_io = true;
  }
  
  // Classify based on connectivity
  if (has_io) {
    return NetType::IO;
  } else {
    return NetType::INTERNAL;
  }
}

double CostEvaluator::computeCost(BStarTree* tree,
                                  double internal_weight,
                                  double io_weight,
                                  double overlap_weight,
                                  double outline_weight)
{
  // Normalize each wirelength class by its baseline so the internal/IO
  // weights trade off comparable magnitudes.  Without this the IO term is
  // typically 1-2% of the total and the weight schedule has no visible
  // effect on the SA landscape.
  const double internal_scale
      = internal_wl_baseline_ > 0.0 ? internal_weight / internal_wl_baseline_
                                    : internal_weight;
  const double io_scale
      = io_wl_baseline_ > 0.0 ? io_weight / io_wl_baseline_ : io_weight;

  double wl_cost = computeWeightedWirelength(tree, internal_scale, io_scale);
  double overlap_cost = computeOverlap(tree);

  // Use placement core dimensions when available (accounts for IO pad ring),
  // otherwise fall back to the full die area.
  const odb::Rect area = placementArea();
  const int max_width = area.dx();
  const int max_height = area.dy();

  double outline_cost = computeOutlinePenalty(tree, max_width, max_height);

  // Normalize penalty areas by the core area so they are dimensionless
  // fractions, matching the normalized wirelength terms.
  const double core_area
      = static_cast<double>(max_width) * static_cast<double>(max_height);
  if (core_area > 0.0) {
    overlap_cost /= core_area;
    outline_cost /= core_area;
  }

  double cost = wl_cost +
                overlap_weight * overlap_cost +
                outline_weight * outline_cost;
  // pack() already keeps macros on the row grid, so only channels can cost.
  if (channel_rule_.min_x > 0 || channel_rule_.min_y > 0) {
    cost += channel_weight_ * computeChannelPenalty(tree);
  }
  // In core half-perimeters per net: at most 1 while each net's soft part
  // lies inside the core.
  if (objective_ == Objective::SPREAD && num_hard_nets_ > 0) {
    const double bound = static_cast<double>(num_hard_nets_)
                         * (static_cast<double>(max_width) + max_height);
    if (bound > 0.0) {
      cost -= spread_weight_ * hard_extension_ / bound;
    }
  }
  return cost;
}

void CostEvaluator::setWirelengthBaselines(double internal_base, double io_base)
{
  internal_wl_baseline_ = internal_base > 0.0 ? internal_base : 0.0;
  io_wl_baseline_ = io_base > 0.0 ? io_base : 0.0;
}

void CostEvaluator::setObjective(Objective objective, double spread_weight)
{
  objective_ = objective;
  spread_weight_ = spread_weight;
}

odb::Rect CostEvaluator::placementArea() const
{
  if (use_placement_core_) {
    return placement_core_;
  }
  return db_->getChip()->getBlock()->getDieArea();
}

void CostEvaluator::collectPlacementBlockages()
{
  placement_blockages_.clear();

  odb::dbBlock* block = db_->getChip()->getBlock();
  const odb::Rect die_area = block->getDieArea();

  for (odb::dbInst* inst : block->getInsts()) {
    if (macro_map_.find(inst) != macro_map_.end()) {
      continue;
    }

    const odb::dbMasterType master_type = inst->getMaster()->getType();
    const bool is_physical_obstacle = inst->isPad() || master_type.isCover()
                                      || master_type == odb::dbMasterType::RING
                                      || inst->isFixed();
    const bool is_placed = inst->isPlaced() || inst->isFixed();
    if (!is_physical_obstacle || !is_placed) {
      continue;
    }

    odb::Rect bbox = inst->getBBox()->getBox();
    if (bbox.intersects(die_area)) {
      placement_blockages_.push_back(bbox);
    }
  }

  for (odb::dbBlockage* blockage : block->getBlockages()) {
    odb::Rect bbox = blockage->getBBox()->getBox();
    if (bbox.intersects(die_area)) {
      placement_blockages_.push_back(bbox);
    }
  }
}

double CostEvaluator::computeWeightedWirelength(BStarTree* tree,
                                                double internal_weight,
                                                double io_weight)
{
  internal_wl_ = 0.0;
  io_wl_ = 0.0;
  internal_soft_wl_ = 0.0;
  io_soft_wl_ = 0.0;
  hard_extension_ = 0.0;

  const odb::Rect core = placementArea();
  for (NetInfo& net_info : nets_) {
    updateNetBoundingBox(net_info, tree, core);
    
    double hpwl = (net_info.max_x - net_info.min_x) + 
                  (net_info.max_y - net_info.min_y);
    
    net_info.hpwl = static_cast<int>(hpwl);
    
    if (net_info.type == NetType::INTERNAL) {
      internal_wl_ += hpwl;
      internal_soft_wl_ += net_info.soft_hpwl;
    } else {
      io_wl_ += hpwl;
      io_soft_wl_ += net_info.soft_hpwl;
    }
    hard_extension_ += net_info.hard_extension;
  }

  if (objective_ == Objective::SPREAD) {
    return internal_weight * internal_soft_wl_ + io_weight * io_soft_wl_;
  }
  return internal_weight * internal_wl_ + io_weight * io_wl_;
}

void CostEvaluator::updateNetBoundingBox(NetInfo& net_info,
                                         BStarTree* tree,
                                         const odb::Rect& core)
{
  net_info.min_x = std::numeric_limits<int>::max();
  net_info.max_x = std::numeric_limits<int>::min();
  net_info.min_y = std::numeric_limits<int>::max();
  net_info.max_y = std::numeric_limits<int>::min();

  const int core_x = use_placement_core_ ? placement_core_.xMin() : 0;
  const int core_y = use_placement_core_ ? placement_core_.yMin() : 0;

  // SPREAD splits the net: soft holds the soft macro centers and fixed
  // points, all adds the hard macro pins.
  const bool spread = objective_ == Objective::SPREAD;
  PointBox soft;
  PointBox all;

  // Movable pins: evaluated from the current B*-tree state, so the SA
  // sees wirelength changes move by move without any DB round-trip.
  const auto& nodes = tree->getNodes();
  for (const auto& [node_id, master_pin] : net_info.macro_pins) {
    if (node_id < 0 || node_id >= static_cast<int>(nodes.size())) {
      continue;
    }
    const BStarNode* node = nodes[node_id].get();
    const odb::Point pos
        = pinPositionFromNode(node, master_pin, core_x, core_y);
    net_info.min_x = std::min(net_info.min_x, pos.x());
    net_info.max_x = std::max(net_info.max_x, pos.x());
    net_info.min_y = std::min(net_info.min_y, pos.y());
    net_info.max_y = std::max(net_info.max_y, pos.y());

    if (!spread) {
      continue;
    }
    if (node->isSoftMacro()) {
      soft.add(pos.x(), pos.y());
      all.add(pos.x(), pos.y());
    } else {
      // Clamped so that leaving the core earns nothing.
      all.add(std::clamp(pos.x(), core.xMin(), core.xMax()),
              std::clamp(pos.y(), core.yMin(), core.yMax()));
    }
  }

  // Fixed connection points (IO pins, pads, placed cells) cached at
  // classification time as a pre-reduced bounding box.
  if (net_info.has_fixed) {
    net_info.min_x = std::min(net_info.min_x, net_info.fixed_min_x);
    net_info.max_x = std::max(net_info.max_x, net_info.fixed_max_x);
    net_info.min_y = std::min(net_info.min_y, net_info.fixed_min_y);
    net_info.max_y = std::max(net_info.max_y, net_info.fixed_max_y);
    if (spread) {
      for (PointBox* box : {&soft, &all}) {
        box->add(net_info.fixed_min_x, net_info.fixed_min_y);
        box->add(net_info.fixed_max_x, net_info.fixed_max_y);
      }
    }
  }

  if (spread) {
    net_info.soft_hpwl = soft.hpwl();
    net_info.hard_extension = all.hpwl() - net_info.soft_hpwl;
  }

  // Handle degenerate case
  if (net_info.min_x == std::numeric_limits<int>::max()) {
    net_info.min_x = 0;
    net_info.max_x = 0;
    net_info.min_y = 0;
    net_info.max_y = 0;
  }
}

double CostEvaluator::computeOverlap(BStarTree* tree)
{
  double total_overlap = 0.0;

  const auto& nodes = tree->getNodes();

  // Tree coordinates are core-relative; blockages are stored in die
  // coordinates.  Shift nodes into the die frame so both comparisons
  // (node-node and node-blockage) are consistent.
  const int core_x = use_placement_core_ ? placement_core_.xMin() : 0;
  const int core_y = use_placement_core_ ? placement_core_.yMin() : 0;

  std::vector<odb::Rect> boxes;
  std::vector<odb::Rect> bodies;
  boxes.reserve(nodes.size());
  bodies.reserve(nodes.size());
  for (const auto& node : nodes) {
    const int x = core_x + node->getX();
    const int y = core_y + node->getY();
    boxes.emplace_back(x, y, x + node->getWidth(), y + node->getHeight());
    const int body_x = core_x + node->getBodyX();
    const int body_y = core_y + node->getBodyY();
    bodies.emplace_back(body_x,
                        body_y,
                        body_x + node->getMacroWidth(),
                        body_y + node->getMacroHeight());
  }

  // Halos may overlap each other, but no body may enter another node's halo
  // box, and no halo a placement blockage.
  for (size_t i = 0; i < nodes.size(); ++i) {
    for (size_t j = i + 1; j < nodes.size(); ++j) {
      const double area = overlapArea(bodies[i], boxes[j])
                          + overlapArea(boxes[i], bodies[j]);
      const bool soft_hard
          = nodes[i]->isSoftMacro() != nodes[j]->isSoftMacro();
      total_overlap += soft_hard ? soft_overlap_weight_ * area : area;
    }

    for (const odb::Rect& blockage : placement_blockages_) {
      total_overlap += overlapArea(boxes[i], blockage);
    }
  }

  return total_overlap;
}

double CostEvaluator::computeOutlinePenalty(BStarTree* tree, 
                                            int max_width, 
                                            int max_height)
{
  double penalty = 0.0;
  for (int subtree = 0; subtree < tree->getNumSubtrees(); ++subtree) {
    int tree_width;
    int tree_height;
    // Add placement core offset to account for IO pad ring if core
    // dimensions are set,
    if (use_placement_core_) {
      tree_width = tree->getSubtreeWidth(subtree) + placement_core_.xMin();
      tree_height = tree->getSubtreeHeight(subtree) + placement_core_.yMin();
    } else {
      tree_width = tree->getSubtreeWidth(subtree);
      tree_height = tree->getSubtreeHeight(subtree);
    }

    if (tree_width > max_width) {
      penalty += static_cast<double>(tree_width - max_width)
                 * static_cast<double>(tree_height);
    }

    if (tree_height > max_height) {
      penalty += static_cast<double>(tree_height - max_height)
                 * static_cast<double>(tree_width);
    }
  }
  
  return penalty;
}

void CostEvaluator::setPlacementCore(const odb::Rect& core)
{
  placement_core_ = core;
  use_placement_core_ = true;
}

void CostEvaluator::setChannelRule(const ChannelRule& rule,
                                   const std::vector<ChannelMacro>& fixed_macros,
                                   double weight)
{
  channel_rule_ = rule;
  fixed_channel_macros_ = fixed_macros;
  channel_weight_ = weight;
}

double CostEvaluator::computeChannelPenalty(const BStarTree* tree) const
{
  return narrowChannelPenalty(
      findNarrowChannels(channelMacros(tree), channel_rule_));
}

std::vector<ChannelMacro> CostEvaluator::channelMacros(
    const BStarTree* tree) const
{
  const int core_x = use_placement_core_ ? placement_core_.xMin() : 0;
  const int core_y = use_placement_core_ ? placement_core_.yMin() : 0;

  std::vector<ChannelMacro> macros;
  macros.reserve(tree->getNodes().size() + fixed_channel_macros_.size());
  for (const auto& node : tree->getNodes()) {
    if (!node->isHardMacro()) {
      continue;
    }
    const Halo& halo = node->getHalo();
    const int x = core_x + node->getX() + halo.left;
    const int y = core_y + node->getY() + halo.bottom;
    ChannelMacro macro;
    macro.body
        = {x, y, x + node->getMacroWidth(), y + node->getMacroHeight()};
    macro.halo = halo;
    macros.push_back(macro);
  }
  macros.insert(
      macros.end(), fixed_channel_macros_.begin(), fixed_channel_macros_.end());
  return macros;
}

}  // namespace pne
