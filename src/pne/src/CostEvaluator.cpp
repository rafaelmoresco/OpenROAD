// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/CostEvaluator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

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
  int max_width;
  int max_height;
  if (use_placement_core_) {
    max_width = placement_core_.dx();
    max_height = placement_core_.dy();
  } else {
    odb::dbBlock* block = db_->getChip()->getBlock();
    odb::Rect die_area = block->getDieArea();
    max_width = die_area.dx();
    max_height = die_area.dy();
  }

  double outline_cost = computeOutlinePenalty(tree, max_width, max_height);

  // Normalize penalty areas by the core area so they are dimensionless
  // fractions, matching the normalized wirelength terms.
  const double core_area
      = static_cast<double>(max_width) * static_cast<double>(max_height);
  if (core_area > 0.0) {
    overlap_cost /= core_area;
    outline_cost /= core_area;
  }

  return wl_cost +
         overlap_weight * overlap_cost +
         outline_weight * outline_cost;
}

void CostEvaluator::setWirelengthBaselines(double internal_base, double io_base)
{
  internal_wl_baseline_ = internal_base > 0.0 ? internal_base : 0.0;
  io_wl_baseline_ = io_base > 0.0 ? io_base : 0.0;
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

double CostEvaluator::computeWirelength(BStarTree* tree)
{
  return computeWeightedWirelength(tree, 1.0, 1.0);
}

double CostEvaluator::computeWeightedWirelength(BStarTree* tree,
                                                double internal_weight,
                                                double io_weight)
{
  internal_wl_ = 0.0;
  io_wl_ = 0.0;
  
  for (NetInfo& net_info : nets_) {
    updateNetBoundingBox(net_info, tree);
    
    double hpwl = (net_info.max_x - net_info.min_x) + 
                  (net_info.max_y - net_info.min_y);
    
    net_info.hpwl = static_cast<int>(hpwl);
    
    if (net_info.type == NetType::INTERNAL) {
      internal_wl_ += hpwl;
    } else {
      io_wl_ += hpwl;
    }
  }
  
  return internal_weight * internal_wl_ + io_weight * io_wl_;
}

void CostEvaluator::updateNetBoundingBox(NetInfo& net_info, BStarTree* tree)
{
  net_info.min_x = std::numeric_limits<int>::max();
  net_info.max_x = std::numeric_limits<int>::min();
  net_info.min_y = std::numeric_limits<int>::max();
  net_info.max_y = std::numeric_limits<int>::min();

  const int core_x = use_placement_core_ ? placement_core_.xMin() : 0;
  const int core_y = use_placement_core_ ? placement_core_.yMin() : 0;

  // Movable pins: evaluated from the current B*-tree state, so the SA
  // sees wirelength changes move by move without any DB round-trip.
  const auto& nodes = tree->getNodes();
  for (const auto& [node_id, master_pin] : net_info.macro_pins) {
    if (node_id < 0 || node_id >= static_cast<int>(nodes.size())) {
      continue;
    }
    const odb::Point pos
        = pinPositionFromNode(nodes[node_id].get(), master_pin, core_x, core_y);
    net_info.min_x = std::min(net_info.min_x, pos.x());
    net_info.max_x = std::max(net_info.max_x, pos.x());
    net_info.min_y = std::min(net_info.min_y, pos.y());
    net_info.max_y = std::max(net_info.max_y, pos.y());
  }

  // Fixed connection points (IO pins, pads, placed cells) cached at
  // classification time as a pre-reduced bounding box.
  if (net_info.has_fixed) {
    net_info.min_x = std::min(net_info.min_x, net_info.fixed_min_x);
    net_info.max_x = std::max(net_info.max_x, net_info.fixed_max_x);
    net_info.min_y = std::min(net_info.min_y, net_info.fixed_min_y);
    net_info.max_y = std::max(net_info.max_y, net_info.fixed_max_y);
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

  // Check pairwise overlaps
  for (size_t i = 0; i < nodes.size(); ++i) {
    int x1 = core_x + nodes[i]->getX();
    int y1 = core_y + nodes[i]->getY();
    int w1 = nodes[i]->getWidth();
    int h1 = nodes[i]->getHeight();

    for (size_t j = i + 1; j < nodes.size(); ++j) {
      int x2 = core_x + nodes[j]->getX();
      int y2 = core_y + nodes[j]->getY();
      int w2 = nodes[j]->getWidth();
      int h2 = nodes[j]->getHeight();
      
      // Compute overlap area
      int overlap_x = std::max(0, std::min(x1 + w1, x2 + w2) - std::max(x1, x2));
      int overlap_y = std::max(0, std::min(y1 + h1, y2 + h2) - std::max(y1, y2));
      
      total_overlap += static_cast<double>(overlap_x)
                       * static_cast<double>(overlap_y);
    }

    for (const odb::Rect& blockage : placement_blockages_) {
      int overlap_x = std::max(0,
                               std::min(x1 + w1, blockage.xMax())
                                   - std::max(x1, blockage.xMin()));
      int overlap_y = std::max(0,
                               std::min(y1 + h1, blockage.yMax())
                                   - std::max(y1, blockage.yMin()));

      total_overlap += static_cast<double>(overlap_x)
                       * static_cast<double>(overlap_y);
    }
  }
  
  return total_overlap;
}

double CostEvaluator::computeOutlinePenalty(BStarTree* tree, 
                                            int max_width, 
                                            int max_height)
{
  int tree_width;
  int tree_height;
  // Add placement core offset to account for IO pad ring if core dimensions are set,
  if (use_placement_core_) {
    tree_width = tree->getWidth() + placement_core_.xMin();
    tree_height = tree->getHeight() + placement_core_.yMin();
  } else {
    tree_width = tree->getWidth();
    tree_height = tree->getHeight();
  }
  
  double penalty = 0.0;
  
  if (tree_width > max_width) {
    penalty += static_cast<double>(tree_width - max_width)
               * static_cast<double>(tree_height);
  }
  
  if (tree_height > max_height) {
    penalty += static_cast<double>(tree_height - max_height)
               * static_cast<double>(tree_width);
  }
  
  return penalty;
}

void CostEvaluator::setPlacementCore(const odb::Rect& core)
{
  placement_core_ = core;
  use_placement_core_ = true;
}

}  // namespace pne
