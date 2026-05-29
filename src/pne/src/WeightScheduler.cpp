// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/WeightScheduler.h"

#include <cmath>
#include <algorithm>

namespace pne {

WeightScheduler::WeightScheduler()
{
  reset();
}

void WeightScheduler::reset()
{
  current_iteration_ = 0;
  current_internal_weight_ = initial_internal_weight_;
  current_io_weight_ = initial_io_weight_;
}

void WeightScheduler::nextIteration()
{
  current_iteration_++;
  
  if (current_iteration_ >= num_iterations_) {
    current_internal_weight_ = final_internal_weight_;
    current_io_weight_ = final_io_weight_;
    return;
  }
  
  double progress = getProgress();
  
  switch (schedule_type_) {
    case ScheduleType::LINEAR:
      current_internal_weight_ = interpolateLinear(
          initial_internal_weight_, final_internal_weight_, progress);
      current_io_weight_ = interpolateLinear(
          initial_io_weight_, final_io_weight_, progress);
      break;
      
    case ScheduleType::EXPONENTIAL:
      current_internal_weight_ = interpolateExponential(
          initial_internal_weight_, final_internal_weight_, progress);
      current_io_weight_ = interpolateExponential(
          initial_io_weight_, final_io_weight_, progress);
      break;
      
    case ScheduleType::STEP:
      current_internal_weight_ = interpolateStep(
          initial_internal_weight_, final_internal_weight_, progress);
      current_io_weight_ = interpolateStep(
          initial_io_weight_, final_io_weight_, progress);
      break;
      
    case ScheduleType::CUSTOM:
      // Can be extended for custom schedules
      current_internal_weight_ = interpolateLinear(
          initial_internal_weight_, final_internal_weight_, progress);
      current_io_weight_ = interpolateLinear(
          initial_io_weight_, final_io_weight_, progress);
      break;
  }
}

double WeightScheduler::getProgress() const
{
  if (num_iterations_ <= 1) {
    return 1.0;
  }
  
  return static_cast<double>(current_iteration_) / 
         static_cast<double>(num_iterations_ - 1);
}

double WeightScheduler::interpolateLinear(double start, double end, double t)
{
  return start + (end - start) * t;
}

double WeightScheduler::interpolateExponential(double start, double end, double t)
{
  // Exponential interpolation: y = start * (end/start)^t
  if (start <= 0 || end <= 0) {
    return interpolateLinear(start, end, t);
  }
  
  return start * std::pow(end / start, t);
}

double WeightScheduler::interpolateStep(double start, double end, double t)
{
  // Step function at midpoint
  if (t < 0.5) {
    return start;
  } else {
    return end;
  }
}

void WeightScheduler::updateWeightsBasedOnIOProportion(double wl_io,
                                                       double wl_internal)
{
  // Compute total wirelength
  double wl_total = wl_io + wl_internal;
  
  if (wl_total <= 0.0) {
    // No nets; use default weights
    current_internal_weight_ = final_internal_weight_;
    current_io_weight_ = final_io_weight_;
    return;
  }
  
  // Compute IO proportion: (WLio * 100) / TWL
  double io_proportion = (wl_io * 100.0) / wl_total;
  
  // Clamp to observed range [kMinIOProportion, kMaxIOProportion]
  if (io_proportion < kMinIOProportion) {
    io_proportion = kMinIOProportion;
  } else if (io_proportion > kMaxIOProportion) {
    io_proportion = kMaxIOProportion;
  }
  
  // Map IO proportion to interpolation parameter [0.0, 1.0]
  // Low proportion (few IO nets) → t ≈ 0 → favor internal weight
  // High proportion (many IO nets) → t ≈ 1 → favor IO weight
  double range = kMaxIOProportion - kMinIOProportion;
  double t = (io_proportion - kMinIOProportion) / range;
  
  // Interpolate weights based on IO proportion
  current_internal_weight_ = interpolateLinear(
      initial_internal_weight_, final_internal_weight_, t);
  current_io_weight_ = interpolateLinear(
      initial_io_weight_, final_io_weight_, t);
}

}  // namespace pne
