// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/SimulatedAnnealing.h"

#include <cmath>
#include <limits>

#include "utl/Logger.h"
#include "pne/CostEvaluator.h"

namespace pne {

SimulatedAnnealing::SimulatedAnnealing(utl::Logger* logger)
    : logger_(logger),
      rng_(std::random_device{}()),
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
  tree->save();
  
  logger_->info(utl::PNE, 31, 
                "Initial cost: {:.2f}, Temperature: {:.2f}",
                current_cost_, current_temp_);
  
  int iteration_in_temp = 0;
  
  while (current_iteration_ < config_.max_iterations &&
         current_temp_ > config_.final_temperature) {
    
    // Perform perturbation
    perturb(tree);
    tree->pack();
    
    double new_cost = cost_function(tree);
    double delta_cost = new_cost - current_cost_;
    
    // Accept or reject
    if (accept(delta_cost)) {
      current_cost_ = new_cost;
      tree->save();
      num_accepted_++;
      
      // Update best solution
      if (new_cost < best_cost_) {
        best_cost_ = new_cost;
        iterations_since_improvement_ = 0;
        
        logger_->info(utl::PNE, 32,
                      "Iteration {}: New best cost {:.2f}",
                      current_iteration_, best_cost_);
      } else {
        iterations_since_improvement_++;
      }
    } else {
      tree->restore();
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
                    "Temp: {:.2f}, Cost: {:.2f}, Accept ratio: {:.2f}%",
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
  tree->restore();
  tree->pack();
  
  logger_->info(utl::PNE, 35,
                "SA completed - Best cost: {:.2f}, Iterations: {}, "
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
  } else if (r < config_.swap_prob + config_.rotate_prob) {
    return PerturbationType::ROTATE;
  } else {
    return PerturbationType::MOVE;
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

}  // namespace pne
