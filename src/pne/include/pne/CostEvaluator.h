// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <vector>
#include <unordered_map>

#include "odb/db.h"
#include "odb/geom.h"
#include "pne/BStarTree.h"

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
  IO,        // Connects macro to IO pad/external pin
  MIXED      // Connects both macros and standard cells
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
};

// Cost evaluation for PineMP placement
class CostEvaluator
{
 public:
  CostEvaluator(odb::dbDatabase* db, 
                sta::dbNetwork* network,
                utl::Logger* logger);
  
  // Net classification
  void classifyNets(const std::vector<odb::dbInst*>& macros);
  
  // Cost computation
  double computeCost(BStarTree* tree,
                     double internal_weight,
                     double io_weight,
                     double overlap_weight,
                     double outline_weight);
  
  // Individual cost components
  double computeWirelength(BStarTree* tree);
  double computeWeightedWirelength(BStarTree* tree,
                                   double internal_weight,
                                   double io_weight);
  double computeOverlap(BStarTree* tree);
  double computeOutlinePenalty(BStarTree* tree, int max_width, int max_height);
  
  // Statistics
  int getNumInternalNets() const { return num_internal_nets_; }
  int getNumIONets() const { return num_io_nets_; }
  int getNumMixedNets() const { return num_mixed_nets_; }
  
  double getInternalWirelength() const { return internal_wl_; }
  double getIOWirelength() const { return io_wl_; }
  double getTotalWirelength() const { return internal_wl_ + io_wl_; }
  
  // Get net information
  const std::vector<NetInfo>& getNets() const { return nets_; }
  
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
  
  // Statistics
  int num_internal_nets_ = 0;
  int num_io_nets_ = 0;
  int num_mixed_nets_ = 0;
  
  double internal_wl_ = 0.0;
  double io_wl_ = 0.0;
  
  // Helper methods
  NetType classifyNet(odb::dbNet* net);
  void collectPlacementBlockages();
  void updateNetBoundingBox(NetInfo& net_info, BStarTree* tree);
  bool isIOPin(odb::dbITerm* iterm);
  bool isIOPin(odb::dbBTerm* bterm);
};

}  // namespace pne
