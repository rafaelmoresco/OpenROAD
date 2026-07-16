// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/SimulatedAnnealing.h"

#include <cmath>
#include <limits>

#include "utl/Logger.h"
#include "pne/CostEvaluator.h"

namespace pne {

SimulatedAnnealing::SimulatedAnnealing(utl::Logger* logger, unsigned int seed)
    : logger_(logger),
      rng_(seed),
      uniform_dist_(0.0, 1.0)
{
}

void SimulatedAnnealing::optimize(BStarTree* tree,
                                  std::function<double(BStarTree*)> cost_function)
{
  logger_->info(utl::PNE, 30, "Starting Simulated Annealing optimization");
  
  // Initialize
  current_temp_ = config_.initial_temperature;
  current_iteration_ = 0;
  num_accepted_ = 0;
  num_rejected_ = 0;
  iterations_since_improvement_ = 0;
  
  // Compute initial cost
  tree->pack();
  current_cost_ = cost_function(tree);
  best_cost_ = current_cost_;
  tree->saveSnapshot(BStarTree::SnapshotSlot::CURRENT);
  tree->saveSnapshot(BStarTree::SnapshotSlot::BEST);

  // Auto-calibrate the initial temperature so that a "typical" uphill move
  // has ~50 % acceptance probability.  With fixed penalty weights the raw
  // cost values can span many orders of magnitude; a temperature that is
  // orders of magnitude below the typical delta_cost turns the SA into a
  // simple greedy descent that cannot escape local optima.
  // The final temperature must be scaled by the same factor: keeping an
  // absolute final temperature above the calibrated initial one would end
  // the anneal before it starts.
  double final_temp = config_.final_temperature;
  if (config_.auto_calibrate_temperature) {
    current_temp_ = calibrateInitialTemperature(tree, cost_function);
    if (config_.initial_temperature > 0.0) {
      final_temp = current_temp_
                   * (config_.final_temperature / config_.initial_temperature);
    }
    logger_->info(utl::PNE, 36,
                  "SA auto-calibrated temperature range: {:.4g} -> {:.4g}",
                  current_temp_, final_temp);
  }

  logger_->info(utl::PNE, 31,
                "Initial cost: {:.4f}, Temperature: {:.4g}",
                current_cost_, current_temp_);

  int iteration_in_temp = 0;

  while (current_iteration_ < config_.max_iterations &&
         current_temp_ > final_temp) {
    
    // Perform perturbation
    perturb(tree);
    tree->pack();
    
    double new_cost = cost_function(tree);
    double delta_cost = new_cost - current_cost_;
    
    // Accept or reject
    if (accept(delta_cost)) {
      current_cost_ = new_cost;
      tree->saveSnapshot(BStarTree::SnapshotSlot::CURRENT);
      num_accepted_++;
      
      // Update best solution
      if (new_cost < best_cost_) {
        best_cost_ = new_cost;
        tree->saveSnapshot(BStarTree::SnapshotSlot::BEST);
        iterations_since_improvement_ = 0;
        
        logger_->info(utl::PNE, 32,
                      "Iteration {}: New best cost {:.4f}",
                      current_iteration_, best_cost_);
      } else {
        iterations_since_improvement_++;
      }
    } else {
      tree->restoreSnapshot(BStarTree::SnapshotSlot::CURRENT);
      num_rejected_++;
      iterations_since_improvement_++;
    }
    
    current_iteration_++;
    iteration_in_temp++;
    
    // Update temperature
    if (iteration_in_temp >= config_.iterations_per_temp) {
      updateTemperature();
      iteration_in_temp = 0;
      
      double accept_ratio = static_cast<double>(num_accepted_) / 
                           (num_accepted_ + num_rejected_);
      
      logger_->info(utl::PNE, 33,
                    "Temp: {:.4g}, Cost: {:.4f}, Accept ratio: {:.2f}%",
                    current_temp_, current_cost_, accept_ratio * 100.0);
    }
    
    // Early stopping
    if (iterations_since_improvement_ >= config_.no_improvement_limit) {
      logger_->info(utl::PNE, 34,
                    "Early stopping: no improvement for {} iterations",
                    config_.no_improvement_limit);
      break;
    }
  }
  
  // Restore best solution
  tree->restoreSnapshot(BStarTree::SnapshotSlot::BEST);
  tree->pack();
  
  logger_->info(utl::PNE, 35,
                "SA completed - Best cost: {:.4f}, Iterations: {}, "
                "Accepted: {}, Rejected: {}",
                best_cost_, current_iteration_, num_accepted_, num_rejected_);
}

void SimulatedAnnealing::perturb(BStarTree* tree)
{
  if (tree->getNumNodes() == 0) {
    return;
  }
  
  PerturbationType type = selectPerturbationType();
  
  switch (type) {
    case PerturbationType::SWAP:
      perturbSwap(tree);
      break;
      
    case PerturbationType::ROTATE:
      perturbRotate(tree);
      break;
      
    case PerturbationType::MOVE:
      perturbMove(tree);
      break;
      
    case PerturbationType::MIXED:
      // Should not reach here as selectPerturbationType handles it
      perturbSwap(tree);
      break;
  }
}

PerturbationType SimulatedAnnealing::selectPerturbationType()
{
  if (config_.perturb_type != PerturbationType::MIXED) {
    return config_.perturb_type;
  }
  
  double r = uniform_dist_(rng_);
  
  if (r < config_.swap_prob) {
    return PerturbationType::SWAP;
  } else if (r < config_.swap_prob + config_.move_prob) {
    return PerturbationType::MOVE;
  } else {
    return PerturbationType::ROTATE;
  }
}

void SimulatedAnnealing::perturbSwap(BStarTree* tree)
{
  int num_nodes = tree->getNumNodes();
  if (num_nodes < 2) {
    return;
  }
  
  int id1 = getRandomNodeId(num_nodes);
  int id2 = getRandomNodeId(num_nodes);
  
  // Ensure different nodes
  while (id2 == id1) {
    id2 = getRandomNodeId(num_nodes);
  }
  
  tree->swapNodes(id1, id2);
}

void SimulatedAnnealing::perturbRotate(BStarTree* tree)
{
  int num_nodes = tree->getNumNodes();
  if (num_nodes == 0) {
    return;
  }
  
  int id = getRandomNodeId(num_nodes);
  tree->rotateNode(id);
}

void SimulatedAnnealing::perturbMove(BStarTree* tree)
{
  int num_nodes = tree->getNumNodes();
  if (num_nodes < 2) {
    return;
  }
  
  int id = getRandomNodeId(num_nodes);
  int new_parent_id = getRandomNodeId(num_nodes);
  
  // Ensure we're not moving to self
  while (new_parent_id == id) {
    new_parent_id = getRandomNodeId(num_nodes);
  }
  
  bool as_left = uniform_dist_(rng_) < 0.5;
  
  tree->moveNode(id, new_parent_id, as_left);
}

int SimulatedAnnealing::getRandomNodeId(int max_id)
{
  std::uniform_int_distribution<int> dist(0, max_id - 1);
  return dist(rng_);
}

bool SimulatedAnnealing::accept(double delta_cost)
{
  // Always accept improvements
  if (delta_cost < 0) {
    return true;
  }
  
  // Accept worse solutions with probability exp(-delta/T)
  double probability = std::exp(-delta_cost / current_temp_);
  double r = uniform_dist_(rng_);
  
  return r < probability;
}

void SimulatedAnnealing::updateTemperature()
{
  current_temp_ *= config_.cooling_rate;
}

double SimulatedAnnealing::calibrateInitialTemperature(
    BStarTree* tree,
    const std::function<double(BStarTree*)>& cost_function)
{
  const int n_samples = config_.calibration_samples;
  double sum_delta = 0.0;
  int n_counted = 0;

  // Record the state we started from so we can restore it afterwards.
  tree->pack();
  double base_cost = cost_function(tree);
  tree->saveSnapshot(BStarTree::SnapshotSlot::CURRENT);

  for (int i = 0; i < n_samples; i++) {
    perturb(tree);
    tree->pack();
    double sample_cost = cost_function(tree);
    sum_delta += std::abs(sample_cost - base_cost);
    n_counted++;

    // Always restore to the exact same starting state so every sample is
    // independent and the tree is not left in a perturbed condition.
    tree->restoreSnapshot(BStarTree::SnapshotSlot::CURRENT);
    tree->pack();
    base_cost = cost_function(tree);
  }

  if (n_counted == 0 || sum_delta == 0.0) {
    return config_.initial_temperature;  // fallback to user-supplied value
  }

  const double avg_delta = sum_delta / n_counted;
  // Set T such that exp(-avg_delta / T) = 0.5, i.e. T = avg_delta / ln(2).
  return avg_delta / std::log(2.0);
}

}  // namespace pne
