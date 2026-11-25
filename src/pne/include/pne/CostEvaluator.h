// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <utility>
#include <vector>
#include <unordered_map>

#include "odb/db.h"
#include "odb/geom.h"
#include "pne/BStarTree.h"
#include "pne/ChannelLegalizer.h"

namespace utl {
class Logger;
}

namespace sta {
class dbNetwork;
}

namespace pne {

// Net classification for weight-based cost computation
enum class NetType {
  INTERNAL,  // Only connects macros
  IO        // Connects macro to IO pad/external pin
};

// Wirelength objective of the annealer.
enum class Objective {
  MIN_WL,  // minimize the HPWL of every net
  SPREAD   // minimize the soft macro part of each net, reward the hard part
};

struct NetInfo {
  odb::dbNet* net;
  NetType type;
  int hpwl;  // Half-perimeter wirelength

  // Bounding box for HPWL computation
  int min_x;
  int max_x;
  int min_y;
  int max_y;

  // Cached geometry so HPWL can be evaluated from B*-tree coordinates
  // during SA, without touching the (stale) DB placement.
  // Movable pins: (tree node id, pin center in master frame, R0).
  std::vector<std::pair<int, odb::Point>> macro_pins;
  // Bounding box of all non-movable connection points (IO pins, pads,
  // placed/fixed cells) in die coordinates, captured at classification
  // time.  Pre-reduced to a box since these never move during SA.
  bool has_fixed = false;
  int fixed_min_x = 0;
  int fixed_max_x = 0;
  int fixed_min_y = 0;
  int fixed_max_y = 0;

  // SPREAD only: HPWL of the soft macro centers and fixed points, and how
  // much the hard macro pins, clamped to the core, add to it.
  int64_t soft_hpwl = 0;
  int64_t hard_extension = 0;
};

// Cost evaluation for PineMP placement
class CostEvaluator
{
 public:
  CostEvaluator(odb::dbDatabase* db, 
                sta::dbNetwork* network,
                utl::Logger* logger);
  
  // Net classification and geometry caching.  Only nets incident to at
  // least one movable macro are kept: all other nets are invariant under
  // macro moves and would only add constant cost and evaluation time.
  // Fixed connection points (IO pins, pads, placed cells) are sampled from
  // the DB at call time, so re-run this after every ppl pin placement.
  void classifyNets(const std::vector<odb::dbInst*>& macros, BStarTree* tree);

  // Cost computation.  Wirelength, overlap and outline components are
  // normalized (wirelength by the baselines below, penalties by the core
  // area) so the weights are comparable across designs and the SA
  // temperature has a consistent scale.
  double computeCost(BStarTree* tree,
                     double internal_weight,
                     double io_weight,
                     double overlap_weight,
                     double outline_weight);

  // Baselines used to normalize the wirelength terms in computeCost.
  // Typically set once from the initial placement so all iterations are
  // scored against the same reference. Non-positive values disable
  // normalization for that term.
  void setWirelengthBaselines(double internal_base, double io_base);
  
  // With SPREAD, computeCost subtracts spread_weight times the hard
  // extension in core half-perimeters per net (README: set_pine_mp_objective).
  void setObjective(Objective objective, double spread_weight);

  // Individual cost components.  The weighted wirelength is that of the
  // minimized part: whole nets, or their soft parts with SPREAD.
  double computeWeightedWirelength(BStarTree* tree,
                                   double internal_weight,
                                   double io_weight);
  // Area of bodies inside other halo boxes (halos may overlap each other) and
  // of halo boxes over placement blockages.  A soft macro over a hard macro
  // counts the soft overlap weight times its area.
  double computeOverlap(BStarTree* tree);
  void setSoftOverlapWeight(double weight) { soft_overlap_weight_ = weight; }
  // Overflow of each subtree's footprint.
  double computeOutlinePenalty(BStarTree* tree, int max_width, int max_height);
  
  // Statistics
  int getNumInternalNets() const { return num_internal_nets_; }
  int getNumIONets() const { return num_io_nets_; }
  
  double getInternalWirelength() const { return internal_wl_; }
  double getIOWirelength() const { return io_wl_; }
  double getTotalWirelength() const { return internal_wl_ + io_wl_; }

  // SPREAD only (DBU): the soft parts of the nets, the hard extension, and
  // the number of nets a hard macro pin can extend.
  double getInternalSoftWirelength() const { return internal_soft_wl_; }
  double getIOSoftWirelength() const { return io_soft_wl_; }
  double getHardExtension() const { return hard_extension_; }
  int getNumHardNets() const { return num_hard_nets_; }
  
  // Get net information
  const std::vector<NetInfo>& getNets() const { return nets_; }
  
  // Set placement core region (area where macros may be placed).
  // When set, the outline penalty uses core dimensions instead of die area.
  void setPlacementCore(const odb::Rect& core);

  // Narrow-channel rule.  Once enabled, computeCost adds weight times the
  // penalty of the tree's hard macros against each other, the fixed macros
  // and the core edge.
  void setChannelRule(const ChannelRule& rule,
                      const std::vector<ChannelMacro>& fixed_macros,
                      double weight);
  double computeChannelPenalty(const BStarTree* tree) const;
  // The tree's hard macros in node order (die coordinates), then the fixed
  // macros.
  std::vector<ChannelMacro> channelMacros(const BStarTree* tree) const;

  // Pads, fixed instances and blockages the macros must not overlap.
  const std::vector<odb::Rect>& getPlacementBlockages() const
  {
    return placement_blockages_;
  }
  
 private:
  odb::dbDatabase* db_;
  sta::dbNetwork* network_;
  utl::Logger* logger_;
  
  // Net information
  std::vector<NetInfo> nets_;
  std::unordered_map<odb::dbNet*, int> net_to_idx_;
  
  // Macro set for classification
  std::unordered_map<odb::dbInst*, bool> macro_map_;
  std::vector<odb::Rect> placement_blockages_;
  
  // Placement core region (area inside IO pad ring)
  odb::Rect placement_core_;
  bool use_placement_core_ = false;

  ChannelRule channel_rule_;
  std::vector<ChannelMacro> fixed_channel_macros_;
  double channel_weight_ = 0.0;
  
  // Statistics
  int num_internal_nets_ = 0;
  int num_io_nets_ = 0;

  double internal_wl_ = 0.0;
  double io_wl_ = 0.0;

  // Normalization baselines for computeCost (raw DBU wirelength).
  double internal_wl_baseline_ = 0.0;
  double io_wl_baseline_ = 0.0;

  Objective objective_ = Objective::MIN_WL;
  double spread_weight_ = 0.0;
  double soft_overlap_weight_ = 1.0;
  // Nets with a hard macro pin and another point, which the hard extension
  // is normalized by.
  int num_hard_nets_ = 0;
  double internal_soft_wl_ = 0.0;
  double io_soft_wl_ = 0.0;
  double hard_extension_ = 0.0;

  // Helper methods
  NetType classifyNet(odb::dbNet* net);
  void collectPlacementBlockages();
  // The placement core, or the die area if none was set.
  odb::Rect placementArea() const;
  // core (die frame) is where SPREAD clamps the hard macro pins.
  void updateNetBoundingBox(NetInfo& net_info,
                            BStarTree* tree,
                            const odb::Rect& core);
};

}  // namespace pne
