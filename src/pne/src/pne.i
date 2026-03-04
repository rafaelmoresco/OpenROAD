// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

%{
#include "ord/OpenRoad.hh"
#include "pne/PineMP.h"
#include "odb/db.h"

namespace ord {
// Defined in OpenRoad.i
pne::PineMP*
getPineMP();
utl::Logger* getLogger();
}

using utl::PNE;
using ord::getPineMP;
%}

%include "../../Exception.i"
%include <std_string.i>

%inline %{

namespace pne {

bool pine_mp_cmd(const int num_threads) {
  auto pine_mp = getPineMP();
  return pine_mp->place(num_threads);
}

void set_debug_cmd(bool debug) {
  auto pine_mp = getPineMP();
  pine_mp->setDebug(debug);
}

void set_num_iterations_cmd(int num_iterations) {
  auto pine_mp = getPineMP();
  pine_mp->setNumIterations(num_iterations);
}

void set_initial_weights_cmd(double internal_weight, double io_weight) {
  auto pine_mp = getPineMP();
  pine_mp->setInitialInternalWeight(internal_weight);
  pine_mp->setInitialIOWeight(io_weight);
}

void set_final_weights_cmd(double internal_weight, double io_weight) {
  auto pine_mp = getPineMP();
  pine_mp->setFinalInternalWeight(internal_weight);
  pine_mp->setFinalIOWeight(io_weight);
}

void set_sa_params_cmd(double initial_temp, double cooling_rate, int max_iterations) {
  auto pine_mp = getPineMP();
  pine_mp->setInitialTemp(initial_temp);
  pine_mp->setCoolingRate(cooling_rate);
  pine_mp->setMaxIterations(max_iterations);
}

void set_pin_strategy_cmd(const char* strategy) {
  auto pine_mp = getPineMP();
  pine_mp->setPinAssignmentStrategy(strategy);
}

} // namespace

%} // inline
