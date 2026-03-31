// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <memory>
#include <vector>

#include "odb/db.h"
#include "utl/Logger.h"

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
  
  // Helper methods
  bool initializePlacement();
  std::vector<odb::dbInst*> collectMacros();
  void buildBStarTree(const std::vector<odb::dbInst*>& macros);
  void runIterativeOptimization();
  void applyFinalPlacement();
  bool runPplIOPlacement(const char* stage_label);
  
  // Get network for timing analysis
  sta::dbNetwork* getNetwork();
};

}  // namespace pne
