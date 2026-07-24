// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/SimulatedAnnealing.h"

#include <algorithm>
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

  // Fast-SA temperature-step state: n counts temperature steps, and the
  // average |cost change| over the current step drives the next temperature.
  int temp_step = 1;
  double step_delta_sum = 0.0;
  int step_delta_count = 0;

  if (config_.use_fast_sa) {
    // Stage 1: T_1 = delta_avg / -ln(P).  P near 1 makes ln(P) a small
    // negative, so T_1 is a high temperature at which almost every uphill
    // move is accepted (random search).
    const double delta_avg = sampleAverageUphillCost(tree, cost_function);
    const double p
        = std::min(0.999999, std::max(1e-6, config_.fast_sa_accept_prob));
    fast_sa_t1_ = (delta_avg > 0.0) ? delta_avg / (-std::log(p))
                                    : config_.initial_temperature;
    current_temp_ = fast_sa_t1_;
    logger_->info(utl::PNE, 37,
                  "Fast-SA schedule: delta_avg={:.4g}, T1={:.4g} "
                  "(P={:.4g}, c={:.4g}, k={})",
                  delta_avg, fast_sa_t1_, config_.fast_sa_accept_prob,
                  config_.fast_sa_c, config_.fast_sa_k);
  } else if (config_.auto_calibrate_temperature) {
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

  // Seed the slack data for the first temperature step's biased moves.
  if (config_.use_slack_moves) {
    computeSlacks(tree);
  }

  int iteration_in_temp = 0;

  // Fast-SA re-heats in stage 3, so it is bounded by the iteration budget and
  // early stopping rather than by a temperature floor.
  while (current_iteration_ < config_.max_iterations
         && (config_.use_fast_sa || current_temp_ > final_temp)) {

    // Perform perturbation
    perturb(tree);
    tree->pack();

    double new_cost = cost_function(tree);
    double delta_cost = new_cost - current_cost_;

    // Track the magnitude of cost changes at this temperature; the Fast-SA
    // schedule uses their average to set the next temperature.
    step_delta_sum += std::abs(delta_cost);
    step_delta_count++;

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
    
    // Advance the temperature at the end of each step.
    if (iteration_in_temp >= config_.iterations_per_temp) {
      if (config_.use_fast_sa) {
        const double delta_cost_avg
            = (step_delta_count > 0) ? step_delta_sum / step_delta_count : 0.0;
        ++temp_step;
        current_temp_ = fastSaTemperature(temp_step, delta_cost_avg);
        step_delta_sum = 0.0;
        step_delta_count = 0;
        // Entering the hill-climbing stage: give it the full no-improvement
        // budget rather than inheriting the greedy stage's exhausted count.
        if (temp_step == config_.fast_sa_k + 1) {
          iterations_since_improvement_ = 0;
        }
      } else {
        updateTemperature();
      }
      iteration_in_temp = 0;

      // Refresh slack for the next step: the arrangement has drifted, so the
      // set of critical blocks and the violating dimension may have changed.
      if (config_.use_slack_moves) {
        computeSlacks(tree);
      }

      double accept_ratio = static_cast<double>(num_accepted_) /
                           (num_accepted_ + num_rejected_);

      logger_->info(utl::PNE, 33,
                    "Temp: {:.4g}, Cost: {:.4f}, Accept ratio: {:.2f}%",
                    current_temp_, current_cost_, accept_ratio * 100.0);
    }

    // Early stopping (disabled when no_improvement_limit <= 0).
    if (config_.no_improvement_limit > 0
        && iterations_since_improvement_ >= config_.no_improvement_limit) {
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

bool SimulatedAnnealing::useSlackBias()
{
  // slack_active_ is set by computeSlacks only while a hard macro overflows
  // the core; when everything fits, moves stay uniform.
  return config_.use_slack_moves && slack_active_
         && uniform_dist_(rng_) < config_.slack_move_prob;
}

void SimulatedAnnealing::perturbSwap(BStarTree* tree)
{
  int num_nodes = tree->getNumNodes();
  if (num_nodes < 2) {
    return;
  }

  int id1;
  int id2;
  if (useSlackBias()) {
    // Move a critical block (sets the violating dimension) into the slot of
    // a block with spare room.
    id1 = selectCriticalNodeId(num_nodes);
    id2 = selectSpaciousNodeId(num_nodes);
  } else {
    id1 = getRandomNodeId(num_nodes);
    id2 = getRandomNodeId(num_nodes);
  }

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

  // Flipping a critical block in the violating dimension is the move most
  // likely to shrink it.
  int id = useSlackBias() ? selectCriticalNodeId(num_nodes)
                          : getRandomNodeId(num_nodes);
  tree->rotateNode(id);
}

void SimulatedAnnealing::perturbMove(BStarTree* tree)
{
  int num_nodes = tree->getNumNodes();
  if (num_nodes < 2) {
    return;
  }
  
  int id;
  int new_parent_id;
  if (useSlackBias()) {
    // Relocate a critical block next to a block that has room to spare.
    id = selectCriticalNodeId(num_nodes);
    new_parent_id = selectSpaciousNodeId(num_nodes);
  } else {
    id = getRandomNodeId(num_nodes);
    new_parent_id = getRandomNodeId(num_nodes);
  }

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

void SimulatedAnnealing::computeSlacks(BStarTree* tree)
{
  const auto& nodes = tree->getNodes();
  const int n = static_cast<int>(nodes.size());
  node_x_slack_.assign(n, 0);
  node_y_slack_.assign(n, 0);
  slack_active_ = false;
  if (n == 0) {
    return;
  }

  // Feasibility gate: slack-biased moves are only useful while a *hard* macro
  // is outside the core.  Soft macros are virtual (no DB blockage, no PDN
  // constraint), and forcing geometric moves once the hard macros already fit
  // only fights wirelength and freezes the search.  Measure the hard-macro
  // bounding box and engage the bias only when it overflows the core.
  const int core_w = tree->getCoreWidth();
  const int core_h = tree->getCoreHeight();
  long long hard_left = std::numeric_limits<long long>::max();
  long long hard_bottom = std::numeric_limits<long long>::max();
  long long hard_right = 0;
  long long hard_top = 0;
  bool has_hard = false;
  for (const auto& node : nodes) {
    if (!node->isHardMacro()) {
      continue;
    }
    has_hard = true;
    hard_left = std::min(hard_left, static_cast<long long>(node->getX()));
    hard_bottom = std::min(hard_bottom, static_cast<long long>(node->getY()));
    hard_right = std::max(hard_right,
                          static_cast<long long>(node->getX()) + node->getWidth());
    hard_top = std::max(hard_top,
                        static_cast<long long>(node->getY()) + node->getHeight());
  }

  if (!has_hard || core_w <= 0 || core_h <= 0) {
    return;  // nothing to keep inside, or no core reference: leave bias off
  }

  const long long hard_w = hard_right - hard_left;
  const long long hard_h = hard_top - hard_bottom;
  const bool x_over = hard_w > core_w;
  const bool y_over = hard_h > core_h;
  if (!x_over && !y_over) {
    return;  // hard macros fit: fall back to uniform moves for wirelength
  }
  slack_active_ = true;
  // Target the hard dimension that is proportionally more over the core.
  slack_target_y_ = (static_cast<double>(hard_h) / core_h)
                    >= (static_cast<double>(hard_w) / core_w);

  // Bounding box of the full packing (hard + soft), taken from the node
  // coordinates (not getWidth()/getHeight(), which can be stale after a
  // rejected move restores a snapshot without re-packing).  Slack itself is
  // computed over every node so soft macros act as real obstacles.
  long long used_right = 0;
  long long used_top = 0;
  for (const auto& node : nodes) {
    used_right = std::max(used_right,
                          static_cast<long long>(node->getX()) + node->getWidth());
    used_top = std::max(used_top,
                        static_cast<long long>(node->getY()) + node->getHeight());
  }

  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) {
    order[i] = i;
  }

  // Horizontal slack: how far each block can shift right before it hits a
  // right-neighbour (a block to its right that overlaps in y) or the used
  // extent.  Process right-to-left so each block's successors are done first
  // (a successor is strictly further right).  slack = latest_left - x.
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return nodes[a]->getX() > nodes[b]->getX();
  });
  std::vector<long long> latest_left(n, 0);
  for (int i : order) {
    const int xi = nodes[i]->getX();
    const int yi = nodes[i]->getY();
    const int wi = nodes[i]->getWidth();
    const int hi = nodes[i]->getHeight();
    long long latest_right = used_right;
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      const int xj = nodes[j]->getX();
      if (xj < xi + wi) {
        continue;  // not to the right of i
      }
      const int yj = nodes[j]->getY();
      const int hj = nodes[j]->getHeight();
      const bool y_overlap = !(yi + hi <= yj || yj + hj <= yi);
      if (y_overlap) {
        latest_right = std::min(latest_right, latest_left[j]);
      }
    }
    latest_left[i] = latest_right - wi;
    node_x_slack_[i] = static_cast<int>(std::max(0LL, latest_left[i] - xi));
  }

  // Vertical slack: symmetric, shifting blocks up.
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return nodes[a]->getY() > nodes[b]->getY();
  });
  std::vector<long long> latest_bottom(n, 0);
  for (int i : order) {
    const int xi = nodes[i]->getX();
    const int yi = nodes[i]->getY();
    const int wi = nodes[i]->getWidth();
    const int hi = nodes[i]->getHeight();
    long long latest_top = used_top;
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      const int yj = nodes[j]->getY();
      if (yj < yi + hi) {
        continue;  // not above i
      }
      const int xj = nodes[j]->getX();
      const int wj = nodes[j]->getWidth();
      const bool x_overlap = !(xi + wi <= xj || xj + wj <= xi);
      if (x_overlap) {
        latest_top = std::min(latest_top, latest_bottom[j]);
      }
    }
    latest_bottom[i] = latest_top - hi;
    node_y_slack_[i] = static_cast<int>(std::max(0LL, latest_bottom[i] - yi));
  }
}

int SimulatedAnnealing::selectCriticalNodeId(int num_nodes)
{
  // Tournament selection: sample a few nodes and keep the least slack in the
  // targeted dimension.  A soft bias keeps the search ergodic.
  const std::vector<int>& slack = slack_target_y_ ? node_y_slack_ : node_x_slack_;
  int best = getRandomNodeId(num_nodes);
  for (int t = 1; t < config_.slack_tournament; ++t) {
    const int cand = getRandomNodeId(num_nodes);
    if (cand < static_cast<int>(slack.size())
        && slack[cand] < slack[best]) {
      best = cand;
    }
  }
  return best;
}

int SimulatedAnnealing::selectSpaciousNodeId(int num_nodes)
{
  const std::vector<int>& slack = slack_target_y_ ? node_y_slack_ : node_x_slack_;
  int best = getRandomNodeId(num_nodes);
  for (int t = 1; t < config_.slack_tournament; ++t) {
    const int cand = getRandomNodeId(num_nodes);
    if (cand < static_cast<int>(slack.size())
        && slack[cand] > slack[best]) {
      best = cand;
    }
  }
  return best;
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

double SimulatedAnnealing::sampleAverageUphillCost(
    BStarTree* tree,
    const std::function<double(BStarTree*)>& cost_function)
{
  const int n_samples = config_.calibration_samples;
  double sum_uphill = 0.0;
  int n_uphill = 0;

  // Record the starting state so every sample is an independent one-move
  // deviation from it.
  tree->pack();
  const double base_cost = cost_function(tree);
  tree->saveSnapshot(BStarTree::SnapshotSlot::CURRENT);

  for (int i = 0; i < n_samples; ++i) {
    perturb(tree);
    tree->pack();
    const double sample_cost = cost_function(tree);
    const double delta = sample_cost - base_cost;
    if (delta > 0.0) {
      sum_uphill += delta;
      ++n_uphill;
    }
    tree->restoreSnapshot(BStarTree::SnapshotSlot::CURRENT);
    tree->pack();
  }

  return (n_uphill > 0) ? sum_uphill / n_uphill : 0.0;
}

double SimulatedAnnealing::fastSaTemperature(int step, double delta_cost) const
{
  // A stalled step (no measurable cost change) would freeze the search;
  // fall back to the stage-1 temperature to keep it moving.
  if (delta_cost <= 0.0) {
    return fast_sa_t1_;
  }
  if (step <= config_.fast_sa_k) {
    // Stage 2: the large c drives the temperature toward zero (pseudo-greedy
    // local search).
    return fast_sa_t1_ * delta_cost / (step * config_.fast_sa_c);
  }
  // Stage 3: dropping c makes the temperature jump back up (re-heat), then
  // decay as 1/n for hill-climbing.
  return fast_sa_t1_ * delta_cost / step;
}

}  // namespace pne
