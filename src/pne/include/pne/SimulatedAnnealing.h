// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <random>
#include <functional>
#include <vector>

#include "pne/BStarTree.h"

namespace utl {
class Logger;
}

namespace pne {

class CostEvaluator;

// Perturbation operators for SA
enum class PerturbationType {
  SWAP,      // Swap two nodes
  ROTATE,    // Flip a macro (MY/MX/R180 group; footprint-preserving)
  MOVE,      // Move a node to different location in tree
  MIXED      // Random selection among all types
};

// SA configuration
struct SAConfig {
  // Temperature schedule
  double initial_temperature = 1000.0;
  double final_temperature = 1.0;
  double cooling_rate = 0.95;
  
  // Iteration limits
  int max_iterations = 10000;
  int iterations_per_temp = 100;
  
  // Perturbation probabilities (used when type is MIXED).
  // Flips (ROTATE) preserve the macro footprint but move its pins, which
  // matters once wirelength responds to orientation.
  double swap_prob = 0.4;
  double rotate_prob = 0.2;
  double move_prob = 0.4;

  // Early stopping
  int no_improvement_limit = 2500;

  // When true the optimizer samples a small number of random moves from the
  // initial state and sets initial_temperature so that a "typical" worsening
  // move has approximately 50 % acceptance probability.  This makes the SA
  // temperature schedule self-consistent with the actual cost magnitudes,
  // which can vary by many orders of magnitude depending on the penalty
  // weights and design size.  The final temperature is scaled by the same
  // factor so the configured initial/final ratio (cooling length) is kept.
  bool auto_calibrate_temperature = false;

  // Number of sample perturbations used during auto-calibration.
  int calibration_samples = 50;

  // --- Fast-SA three-stage schedule (Chen & Chang, ISPD'05 / TCAD'06) ---
  // When enabled, the geometric cooling above is replaced by the adaptive
  // three-stage schedule.  With n the temperature-step index (each step is a
  // batch of iterations_per_temp moves) and delta_cost the average |cost
  // change| observed over the previous step:
  //   n = 1        : T_1 = delta_avg / -ln(P)      (high-temp random search)
  //   2 <= n <= k  : T_n = T_1 * delta_cost / (n*c)  (pseudo-greedy dive)
  //   n > k        : T_n = T_1 * delta_cost / n      (re-heat + hill-climb)
  // The large c drives stage 2 toward a greedy local search; removing it at
  // stage 3 makes the temperature jump back up to escape the local minimum,
  // then decay as 1/n.  delta_cost tracks the current landscape, so the
  // schedule self-adapts instead of following a fixed ratio.
  bool use_fast_sa = false;
  double fast_sa_accept_prob = 0.99;  // P: stage-1 uphill acceptance for T_1
  double fast_sa_c = 100.0;           // c: stage-2 temperature suppression
  int fast_sa_k = 7;                  // k: stage-2 -> stage-3 boundary

  // --- Slack-based move selection (Adya & Markov, ParquetFP) ---
  // Each block has a "slack" in x and y: how far it can shift toward the
  // right / top before it hits a neighbour or the used extent.  Zero-slack
  // blocks are critical -- they set the floorplan's width / height.  A
  // fraction of moves are biased to pick a critical block in whichever
  // dimension is currently most over the core outline and relocate it toward
  // a block with spare room, which packs the violating dimension faster than
  // uniform-random moves.  The remaining moves stay uniform for ergodicity.
  bool use_slack_moves = true;
  double slack_move_prob = 0.5;  // fraction of moves that are slack-biased
  int slack_tournament = 3;      // tournament size for the soft slack bias

  // Perturbation type
  PerturbationType perturb_type = PerturbationType::MIXED;
};

// Simulated Annealing optimizer
class SimulatedAnnealing
{
 public:
  SimulatedAnnealing(utl::Logger* logger, unsigned int seed = 42);
  
  // Configuration
  void setConfig(const SAConfig& config) { config_ = config; }
  SAConfig getConfig() const { return config_; }
  
  // Main optimization function
  // cost_function: function that computes cost given the tree
  void optimize(BStarTree* tree,
                std::function<double(BStarTree*)> cost_function);
  
  // Get statistics
  double getBestCost() const { return best_cost_; }
  double getCurrentCost() const { return current_cost_; }
  double getCurrentTemperature() const { return current_temp_; }
  int getCurrentIteration() const { return current_iteration_; }
  int getNumAccepted() const { return num_accepted_; }
  int getNumRejected() const { return num_rejected_; }
  
 private:
  utl::Logger* logger_;
  SAConfig config_;
  
  // State
  double best_cost_ = 0.0;
  double current_cost_ = 0.0;
  double current_temp_ = 0.0;
  int current_iteration_ = 0;
  int num_accepted_ = 0;
  int num_rejected_ = 0;
  int iterations_since_improvement_ = 0;

  // Fast-SA stage-1 temperature (T_1), reused as the scale for stages 2-3.
  double fast_sa_t1_ = 0.0;

  // Slack-based move state (per-node slack, indexed by node id, refreshed at
  // each temperature step).  slack_target_y_ selects which dimension the
  // biased moves try to shrink.
  std::vector<int> node_x_slack_;
  std::vector<int> node_y_slack_;
  bool slack_target_y_ = false;

  // Random number generation
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_dist_;

  // Helper methods
  void perturb(BStarTree* tree);
  bool accept(double delta_cost);
  void updateTemperature();
  double calibrateInitialTemperature(
      BStarTree* tree,
      const std::function<double(BStarTree*)>& cost_function);
  // Average uphill cost of a batch of random moves from the current state;
  // seeds the Fast-SA stage-1 temperature T_1.
  double sampleAverageUphillCost(
      BStarTree* tree,
      const std::function<double(BStarTree*)>& cost_function);
  // Fast-SA temperature for step index n given the average |cost change|
  // measured over the previous step.
  double fastSaTemperature(int step, double delta_cost) const;

  // Refresh per-node x/y slack from the current packing and pick the
  // dimension (x or y) currently most over the core outline to shrink.
  void computeSlacks(BStarTree* tree);
  // Tournament pick of a node with low slack (critical) or high slack
  // (spacious) in the currently targeted dimension.
  int selectCriticalNodeId(int num_nodes);
  int selectSpaciousNodeId(int num_nodes);
  // True with probability slack_move_prob when slack data is available.
  bool useSlackBias();

  PerturbationType selectPerturbationType();
  void perturbSwap(BStarTree* tree);
  void perturbRotate(BStarTree* tree);
  void perturbMove(BStarTree* tree);
  
  int getRandomNodeId(int max_id);
};

}  // namespace pne
