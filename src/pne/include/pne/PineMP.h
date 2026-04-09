// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "odb/db.h"
#include "utl/Logger.h"
#include "pne/BStarTree.h"

namespace odb {
class dbDatabase;
class dbInst;
class dbOrientType;
}  // namespace odb

namespace sta {
class dbNetwork;
class dbSta;
}  // namespace sta

namespace par {
class PartitionMgr;
}

namespace ppl {
class IOPlacer;
}

namespace pne {

class BStarTree;
class CostEvaluator;
class WeightScheduler;
class PinAssigner;
class SimulatedAnnealing;
struct SAConfig;

// PineMP main class - orchestrates iterative co-optimization
class PineMP
{
public:
  PineMP(odb::dbDatabase* db, utl::Logger* logger, ppl::IOPlacer* io_placer);
  ~PineMP();

  // Main placement function
  bool place(int num_threads = 1);
  
  // Set STA engine for timing-aware features (optional)
  void setSta(sta::dbSta* sta) { sta_ = sta; }

  // Configuration
  void setDebug(bool debug) { debug_ = debug; }
  bool getDebug() const { return debug_; }
  
  // Weight schedule configuration
  void setNumIterations(int num) { num_iterations_ = num; }
  void setInitialInternalWeight(double weight) { initial_internal_weight_ = weight; }
  void setInitialIOWeight(double weight) { initial_io_weight_ = weight; }
  void setFinalInternalWeight(double weight) { final_internal_weight_ = weight; }
  void setFinalIOWeight(double weight) { final_io_weight_ = weight; }
  
  // SA configuration
  void setInitialTemp(double temp) { initial_temp_ = temp; }
  void setCoolingRate(double rate) { cooling_rate_ = rate; }
  void setMaxIterations(int max_iter) { max_sa_iterations_ = max_iter; }
  
  // Pin assignment strategy
  void setPinAssignmentStrategy(const std::string& strategy);

  // Halo configuration.  Horizontal halo is applied to left/right sides,
  // vertical halo to top/bottom sides.
  void setHalo(int halo_x, int halo_y) { halo_x_ = halo_x; halo_y_ = halo_y; }
  void setPinAwareHalo(bool enable) { pin_aware_halo_ = enable; }

  // Per-macro halo override (4-sided: left, bottom, right, top in DBU).
  void setMacroHalo(const std::string& macro_name, const Halo& halo);

private:
  utl::Logger* logger_;
  odb::dbDatabase* db_;
  sta::dbSta* sta_ = nullptr;
  
  bool debug_ = false;
  
  // Core components
  std::unique_ptr<BStarTree> tree_;
  std::unique_ptr<CostEvaluator> cost_evaluator_;
  std::unique_ptr<WeightScheduler> weight_scheduler_;
  std::unique_ptr<PinAssigner> pin_assigner_;
  std::unique_ptr<SimulatedAnnealing> sa_optimizer_;
  ppl::IOPlacer* io_placer_ = nullptr;
  
  // Configuration
  int num_iterations_ = 5;
  double initial_internal_weight_ = 0.8;
  double initial_io_weight_ = 0.2;
  double final_internal_weight_ = 0.5;
  double final_io_weight_ = 0.5;
  
  // SA parameters
  double initial_temp_ = 1000.0;
  double cooling_rate_ = 0.95;
  int max_sa_iterations_ = 10000;
  
  // Halo parameters (DBU)
  int halo_x_ = 0;
  int halo_y_ = 0;
  bool pin_aware_halo_ = true;  // Default: only add halo on sides with pins

  // Per-macro halo overrides stored until tree is built
  std::unordered_map<std::string, Halo> macro_halo_overrides_;
  
  // Placement region (core area inside IO pad ring)
  odb::Rect placement_core_;
  
  // Helper methods
  bool initializePlacement();
  void computePlacementRegion();
  void applyPlacementWithOffset();
  std::vector<odb::dbInst*> collectMacros();
  void buildBStarTree(const std::vector<odb::dbInst*>& macros);
  void applyHalos();
  void runIterativeOptimization();
  void applyFinalPlacement();
  bool runPplIOPlacement(const char* stage_label);
  
  // Get network for timing analysis
  sta::dbNetwork* getNetwork();
};

}  // namespace pne
