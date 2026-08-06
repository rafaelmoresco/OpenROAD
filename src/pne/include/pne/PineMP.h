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
#include "pne/RecursivePartitioner.h"
#include "pne/SoftMacro.h"

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
class SimulatedAnnealing;
struct SAConfig;

// PineMP main class - orchestrates iterative co-optimization
class PineMP
{
public:
  PineMP(odb::dbDatabase* db,
         utl::Logger* logger,
         ppl::IOPlacer* io_placer,
         par::PartitionMgr* partition_mgr = nullptr);
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
  
  // Soft macro configuration.
  // Call set_pine_mp_soft_macros then pine_mp to co-place stdcell
  // clusters (from a prior triton_part run) with hard macros.
  void enableSoftMacros(bool enable) { use_soft_macros_ = enable; }
  bool softMacrosEnabled() const { return use_soft_macros_; }
  void setSoftMacroUtilization(double u) { soft_macro_utilization_ = u; }
  void setSoftMacroAspectRatio(double r) { soft_macro_aspect_ratio_ = r; }
  // Report statistics for the current set of soft macros (if any).
  void reportSoftMacros() const;

  // Internal recursive-bisection partitioning (replaces the external
  // triton_part call in the flow).  Enabled by default when soft macros
  // are used; disable to consume pre-existing partition_id properties.
  void enableInternalPartitioning(bool enable)
  {
    use_internal_partitioning_ = enable;
  }
  // Split until at least this many partitions exist (0 = no count target).
  void setPartitionTarget(int num) { partition_target_ = num; }
  // Keep splitting partitions whose soft-macro footprint would exceed this
  // fraction of the core area (0 = no size ceiling).
  void setPartitionMaxAreaFraction(double f)
  {
    partition_max_area_fraction_ = f;
  }
  void setPartitionMinCells(int num) { partition_min_cells_ = num; }
  void setPartitionSeed(int seed) { partition_seed_ = seed; }
  // Report the bisection hierarchy of the last pine_mp run.
  void reportPartitionTree() const;

  // Halo configuration.  Horizontal halo is applied to left/right sides,
  // vertical halo to top/bottom sides.
  void setHalo(int halo_x, int halo_y) { halo_x_ = halo_x; halo_y_ = halo_y; }
  void setPinAwareHalo(bool enable) { pin_aware_halo_ = enable; }

  // Per-macro halo override (4-sided: left, bottom, right, top in DBU).
  void setMacroHalo(const std::string& macro_name, const Halo& halo);

  // Adaptive weight scheduling based on IO wirelength proportion.
  // When enabled, weights are dynamically adjusted based on the ratio of
  // IO wirelength to total wirelength, rather than using a fixed schedule.
  void enableAdaptiveIOWeighting(bool enable) { use_adaptive_io_weighting_ = enable; }
  bool adaptiveIOWeightingEnabled() const { return use_adaptive_io_weighting_; }

  // Four-corner anchoring: after each SA run, re-evaluate the packed layout
  // compacted toward each of the four core corners and keep the lowest-cost
  // one.  Counteracts the bottom-left bias inherent to B*-tree packing by
  // aligning the macro cluster with the fixed IO pins.
  void enableCornerAnchoring(bool enable) { use_corner_anchoring_ = enable; }
  bool cornerAnchoringEnabled() const { return use_corner_anchoring_; }

  // Fast-SA three-stage annealing schedule (Chen & Chang).  When enabled,
  // each SA run uses the adaptive high-temp / pseudo-greedy / re-heat
  // schedule instead of geometric cooling with 50%-acceptance calibration.
  void enableFastSA(bool enable) { use_fast_sa_ = enable; }
  bool fastSAEnabled() const { return use_fast_sa_; }
  void setFastSAParams(double accept_prob, double c, int k)
  {
    fast_sa_accept_prob_ = accept_prob;
    fast_sa_c_ = c;
    fast_sa_k_ = k;
  }

  // Slack-based move selection (Adya & Markov): bias a fraction of SA moves
  // toward blocks critical in the dimension over the outline.
  void enableSlackMoves(bool enable) { use_slack_moves_ = enable; }
  bool slackMovesEnabled() const { return use_slack_moves_; }
  void setSlackMoveProbability(double prob) { slack_move_prob_ = prob; }

private:
  utl::Logger* logger_;
  odb::dbDatabase* db_;
  sta::dbSta* sta_ = nullptr;
  
  bool debug_ = false;
  
  // Core components
  std::unique_ptr<BStarTree> tree_;
  std::unique_ptr<CostEvaluator> cost_evaluator_;
  std::unique_ptr<WeightScheduler> weight_scheduler_;
  std::unique_ptr<SimulatedAnnealing> sa_optimizer_;
  ppl::IOPlacer* io_placer_ = nullptr;

  // Soft macro support
  bool use_soft_macros_ = false;
  double soft_macro_utilization_ = 0.7;
  double soft_macro_aspect_ratio_ = 1.0;
  std::unique_ptr<SoftMacroMgr> soft_macro_mgr_;

  // Internal partitioning (recursive bisection via the par module)
  par::PartitionMgr* partition_mgr_ = nullptr;
  bool use_internal_partitioning_ = true;
  int partition_target_ = 0;              // 0 = no count target
  double partition_max_area_fraction_ = 0.0;  // 0 = no size ceiling
  int partition_min_cells_ = 50;
  int partition_seed_ = 1;
  // Bisection hierarchy of the last run, kept for reporting and future
  // hierarchical placement.
  std::vector<PartitionTreeNode> partition_tree_;
  
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

  // Adaptive weight scheduling
  bool use_adaptive_io_weighting_ = false;  // Use IO proportion to drive weights

  // Four-corner anchoring
  bool use_corner_anchoring_ = true;

  // Fast-SA schedule.  DISABLED by default: the literal Chen & Chang
  // schedule is mis-scaled for PineMP's penalty-heavy normalized cost.
  // T1 = delta_avg / -ln(0.99) ~ 100*delta_avg is very hot, and the stage-3
  // 1/n decay is far too slow to reach the low temperatures needed to pack
  // tightly (the geometric schedule cools ~1000x below delta_avg; Fast-SA's
  // 1/n barely reaches delta_avg within the iteration budget).  The stiff
  // overlap/outline penalties also keep the running delta_cost large, so the
  // self-adaptive schedule never cools -- the search random-walks at high
  // temperature and the macro cluster overflows the core outline.  Kept
  // behind the flag for research; the geometric auto-calibrated schedule is
  // the working default.
  bool use_fast_sa_ = false;
  double fast_sa_accept_prob_ = 0.99;
  double fast_sa_c_ = 100.0;
  int fast_sa_k_ = 7;

  // Slack-based move selection
  bool use_slack_moves_ = true;
  double slack_move_prob_ = 0.5;

  // Per-macro halo overrides stored until tree is built
  std::unordered_map<std::string, Halo> macro_halo_overrides_;
  
  // Placement region (core area inside IO pad ring)
  odb::Rect placement_core_;
  
  // Parallel multi-start: number of independently seeded SA chains per
  // outer iteration (from pine_mp -num_threads); the best chain wins.
  int num_threads_ = 1;

  // Helper methods
  bool initializePlacement();
  void computePlacementRegion();
  void applyPlacementWithOffset();
  void enforceBoundsCompliance();
  std::vector<odb::dbInst*> collectMacros();
  void buildBStarTree(const std::vector<odb::dbInst*>& macros);
  void attachSoftMacros();
  void applyHalos();
  void runIterativeOptimization();
  void applyFinalPlacement();
  bool runPplIOPlacement(const char* stage_label);
  // Try all four corner anchorings on the current tree and commit the one
  // with the lowest cost under the given (validation) objective.
  void selectBestAnchor(double internal_weight,
                        double io_weight,
                        double overlap_weight,
                        double outline_weight);
  // Run the SA for one outer iteration under the given weights.  With
  // num_threads_ > 1, runs that many independently seeded chains in
  // parallel on cloned trees and adopts the lowest-cost result; with 1,
  // runs the legacy single persistent chain (bit-compatible with previous
  // behavior).
  void runMultiStartSA(double internal_weight,
                       double io_weight,
                       double overlap_weight,
                       double outline_weight,
                       int iteration);
  
  // Get network for timing analysis
  sta::dbNetwork* getNetwork();
};

}  // namespace pne
