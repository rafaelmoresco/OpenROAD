// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <vector>
#include <unordered_map>

#include "odb/db.h"

namespace utl {
class Logger;
}

namespace pne {

// Pin assignment strategy
enum class PinAssignmentStrategy {
  UNIFORM,        // Uniform distribution along boundary
  CONNECTIVITY,   // Connectivity-based clustering
  HUNGARAIN,     // Hungarian algorithm for optimal assignment
  RANDOM          // Random assignment
};

// Pin location on macro boundary
struct PinLocation {
  int x;
  int y;
  odb::dbPlacementStatus status;
};

// Pin assignment engine for macro pins
class PinAssigner
{
 public:
  PinAssigner(odb::dbDatabase* db, utl::Logger* logger);
  
  // Configuration
  void setStrategy(PinAssignmentStrategy strategy) { strategy_ = strategy; }
  PinAssignmentStrategy getStrategy() const { return strategy_; }
  
  // Main assignment function
  void assignPins(const std::vector<odb::dbInst*>& macros);
  
  // Assign pins for a single macro
  void assignMacroPins(odb::dbInst* macro);
  
  // Reset to default pin positions
  void resetPins(const std::vector<odb::dbInst*>& macros);
  
 private:
  odb::dbDatabase* db_;
  utl::Logger* logger_;
  
  PinAssignmentStrategy strategy_ = PinAssignmentStrategy::UNIFORM;
  
  // Helper methods
  void assignUniform(odb::dbInst* macro);
  void assignConnectivity(odb::dbInst* macro);
  void assignRandom(odb::dbInst* macro);
  
  // Get pins for a macro
  std::vector<odb::dbITerm*> getMacroPins(odb::dbInst* macro);
  
  // Compute boundary positions
  enum class Side { LEFT, RIGHT, TOP, BOTTOM };
  PinLocation computePinLocation(odb::dbInst* macro, 
                                Side side, 
                                double fraction);
  
  // Connectivity analysis
  struct PinConnectivity {
    odb::dbITerm* iterm;
    double center_x;  // Weighted center of connected pins
    double center_y;
    int num_connections;
  };
  
  std::vector<PinConnectivity> analyzeConnectivity(odb::dbInst* macro);
  
  // Assign pin to boundary based on connectivity
  Side selectBoundarySide(odb::dbInst* macro, 
                         double target_x, 
                         double target_y);
};

}  // namespace pne
