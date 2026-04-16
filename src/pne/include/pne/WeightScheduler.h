// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

namespace pne {

// Weight scheduler for adaptive weight tuning during iterative co-optimization
class WeightScheduler
{
 public:
  WeightScheduler();
  
  // Configuration
  void setInitialInternalWeight(double weight) { initial_internal_weight_ = weight; }
  void setInitialIOWeight(double weight) { initial_io_weight_ = weight; }
  void setFinalInternalWeight(double weight) { final_internal_weight_ = weight; }
  void setFinalIOWeight(double weight) { final_io_weight_ = weight; }
  void setNumIterations(int num) { num_iterations_ = num; }
  
  // Weight schedule types
  enum class ScheduleType {
    LINEAR,      // Linear interpolation
    EXPONENTIAL, // Exponential decay/growth
    STEP,        // Step function
    CUSTOM       // User-defined schedule
  };
  
  void setScheduleType(ScheduleType type) { schedule_type_ = type; }
  
  // Initialize scheduler
  void reset();
  
  // Get current weights
  double getInternalWeight() const { return current_internal_weight_; }
  double getIOWeight() const { return current_io_weight_; }
  double getOverlapWeight() const { return overlap_weight_; }
  double getOutlineWeight() const { return outline_weight_; }
  
  // Update weights for next iteration
  void nextIteration();
  
  // Query state
  int getCurrentIteration() const { return current_iteration_; }
  int getNumIterations() const { return num_iterations_; }
  bool isDone() const { return current_iteration_ >= num_iterations_; }
  
  // Get progress ratio [0.0, 1.0]
  double getProgress() const;
  
  // Fixed weights for overlap and outline penalty
  void setOverlapWeight(double weight) { overlap_weight_ = weight; }
  void setOutlineWeight(double weight) { outline_weight_ = weight; }
  
 private:
  // Configuration
  double initial_internal_weight_ = 0.8;
  double initial_io_weight_ = 0.2;
  double final_internal_weight_ = 0.5;
  double final_io_weight_ = 0.5;
  int num_iterations_ = 5;
  
  ScheduleType schedule_type_ = ScheduleType::LINEAR;
  
  // Current state
  int current_iteration_ = 0;
  double current_internal_weight_ = 0.8;
  double current_io_weight_ = 0.2;
  
  // Fixed weights
  double overlap_weight_ = 1e7;  // High penalty for overlaps
  double outline_weight_ = 1e7;  // High penalty for outline violations
  
  // Helper methods
  double interpolateLinear(double start, double end, double t);
  double interpolateExponential(double start, double end, double t);
  double interpolateStep(double start, double end, double t);
};

}  // namespace pne
