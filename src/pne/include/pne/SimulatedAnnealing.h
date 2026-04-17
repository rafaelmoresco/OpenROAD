// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <random>
#include <functional>

#include "pne/BStarTree.h"

namespace utl {
class Logger;
}

namespace pne {

class CostEvaluator;

// Perturbation operators for SA
enum class PerturbationType {
  SWAP,      // Swap two nodes
  ROTATE,    // Rotate a macro
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
  
  // Perturbation probabilities (used when type is MIXED)
  double swap_prob = 0.5;
  double rotate_prob = 0;
  double move_prob = 0.5;
  
  // Early stopping
  int no_improvement_limit = 1000;
  
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
  
  // Random number generation
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_dist_;
  
  // Helper methods
  void perturb(BStarTree* tree);
  bool accept(double delta_cost);
  void updateTemperature();
  
  PerturbationType selectPerturbationType();
  void perturbSwap(BStarTree* tree);
  void perturbRotate(BStarTree* tree);
  void perturbMove(BStarTree* tree);
  
  int getRandomNodeId(int max_id);
};

}  // namespace pne
