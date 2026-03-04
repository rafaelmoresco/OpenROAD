// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/CostEvaluator.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "utl/Logger.h"
#include "sta/Network.hh"

namespace pne {

CostEvaluator::CostEvaluator(odb::dbDatabase* db,
                             sta::dbNetwork* network,
                             utl::Logger* logger)
    : db_(db), network_(network), logger_(logger)
{
}

void CostEvaluator::classifyNets(const std::vector<odb::dbInst*>& macros)
{
  // Build macro lookup map
  macro_map_.clear();
  for (odb::dbInst* macro : macros) {
    macro_map_[macro] = true;
  }
  
  // Get block and iterate through all nets
  odb::dbBlock* block = db_->getChip()->getBlock();
  nets_.clear();
  net_to_idx_.clear();
  
  num_internal_nets_ = 0;
  num_io_nets_ = 0;
  num_mixed_nets_ = 0;
  
  for (odb::dbNet* net : block->getNets()) {
    // Skip special nets (power, ground, clock)
    if (net->getSigType() == odb::dbSigType::POWER ||
        net->getSigType() == odb::dbSigType::GROUND ||
        net->getSigType() == odb::dbSigType::CLOCK) {
      continue;
    }
    
    NetType type = classifyNet(net);
    
    NetInfo info;
    info.net = net;
    info.type = type;
    info.hpwl = 0;
    info.min_x = std::numeric_limits<int>::max();
    info.max_x = std::numeric_limits<int>::min();
    info.min_y = std::numeric_limits<int>::max();
    info.max_y = std::numeric_limits<int>::min();
    
    net_to_idx_[net] = nets_.size();
    nets_.push_back(info);
    
    switch (type) {
      case NetType::INTERNAL:
        num_internal_nets_++;
        break;
      case NetType::IO:
        num_io_nets_++;
        break;
      case NetType::MIXED:
        num_mixed_nets_++;
        break;
    }
  }
  
  logger_->info(utl::PNE, 10, 
                "Net classification - Internal: {}, IO: {}, Mixed: {}",
                num_internal_nets_, num_io_nets_, num_mixed_nets_);
}

NetType CostEvaluator::classifyNet(odb::dbNet* net)
{
  bool has_macro = false;
  bool has_io = false;
  bool has_std_cell = false;
  
  // Check instance terminals
  for (odb::dbITerm* iterm : net->getITerms()) {
    odb::dbInst* inst = iterm->getInst();
    
    if (macro_map_.find(inst) != macro_map_.end()) {
      has_macro = true;
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
  } else if (has_macro && !has_std_cell) {
    return NetType::INTERNAL;
  } else {
    return NetType::MIXED;
  }
}

double CostEvaluator::computeCost(BStarTree* tree,
                                  double internal_weight,
                                  double io_weight,
                                  double overlap_weight,
                                  double outline_weight)
{
  double wl_cost = computeWeightedWirelength(tree, internal_weight, io_weight);
  double overlap_cost = computeOverlap(tree);
  
  // Get die area for outline penalty
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::Rect die_area = block->getDieArea();
  int max_width = die_area.dx();
  int max_height = die_area.dy();
  
  double outline_cost = computeOutlinePenalty(tree, max_width, max_height);
  
  return wl_cost + 
         overlap_weight * overlap_cost + 
         outline_weight * outline_cost;
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
    } else if (net_info.type == NetType::IO) {
      io_wl_ += hpwl;
    } else {
      // Split mixed nets
      internal_wl_ += hpwl * 0.5;
      io_wl_ += hpwl * 0.5;
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
  
  odb::dbNet* net = net_info.net;
  
  // Process instance terminals
  for (odb::dbITerm* iterm : net->getITerms()) {
    odb::dbInst* inst = iterm->getInst();
    
    // Get pin location
    int x, y;
    if (iterm->getAvgXY(&x, &y)) {
      net_info.min_x = std::min(net_info.min_x, x);
      net_info.max_x = std::max(net_info.max_x, x);
      net_info.min_y = std::min(net_info.min_y, y);
      net_info.max_y = std::max(net_info.max_y, y);
    } else {
      // Use instance location as approximation
      int inst_x, inst_y;
      inst->getLocation(inst_x, inst_y);
      net_info.min_x = std::min(net_info.min_x, inst_x);
      net_info.max_x = std::max(net_info.max_x, inst_x);
      net_info.min_y = std::min(net_info.min_y, inst_y);
      net_info.max_y = std::max(net_info.max_y, inst_y);
    }
  }
  
  // Process block terminals (IO pins)
  for (odb::dbBTerm* bterm : net->getBTerms()) {
    for (odb::dbBPin* bpin : bterm->getBPins()) {
      for (odb::dbBox* box : bpin->getBoxes()) {
        odb::Rect rect = box->getBox();
        int x = (rect.xMin() + rect.xMax()) / 2;
        int y = (rect.yMin() + rect.yMax()) / 2;
        
        net_info.min_x = std::min(net_info.min_x, x);
        net_info.max_x = std::max(net_info.max_x, x);
        net_info.min_y = std::min(net_info.min_y, y);
        net_info.max_y = std::max(net_info.max_y, y);
      }
    }
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
  
  // Check pairwise overlaps
  for (size_t i = 0; i < nodes.size(); ++i) {
    int x1 = nodes[i]->getX();
    int y1 = nodes[i]->getY();
    int w1 = nodes[i]->getWidth();
    int h1 = nodes[i]->getHeight();
    
    for (size_t j = i + 1; j < nodes.size(); ++j) {
      int x2 = nodes[j]->getX();
      int y2 = nodes[j]->getY();
      int w2 = nodes[j]->getWidth();
      int h2 = nodes[j]->getHeight();
      
      // Compute overlap area
      int overlap_x = std::max(0, std::min(x1 + w1, x2 + w2) - std::max(x1, x2));
      int overlap_y = std::max(0, std::min(y1 + h1, y2 + h2) - std::max(y1, y2));
      
      total_overlap += overlap_x * overlap_y;
    }
  }
  
  return total_overlap;
}

double CostEvaluator::computeOutlinePenalty(BStarTree* tree, 
                                            int max_width, 
                                            int max_height)
{
  int tree_width = tree->getWidth();
  int tree_height = tree->getHeight();
  
  double penalty = 0.0;
  
  if (tree_width > max_width) {
    penalty += (tree_width - max_width) * tree_height;
  }
  
  if (tree_height > max_height) {
    penalty += (tree_height - max_height) * tree_width;
  }
  
  return penalty;
}

bool CostEvaluator::isIOPin(odb::dbITerm* iterm)
{
  // This is an instance terminal, not an IO pin
  return false;
}

bool CostEvaluator::isIOPin(odb::dbBTerm* bterm)
{
  // Block terminals are IO pins
  return true;
}

}  // namespace pne
