// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

%{
#include "ord/OpenRoad.hh"
#include "pne/PineMP.h"
#include "pne/BStarTree.h"
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

void set_halo_cmd(int halo_x, int halo_y, bool pin_aware) {
  auto pine_mp = getPineMP();
  pine_mp->setHalo(halo_x, halo_y);
  pine_mp->setPinAwareHalo(pin_aware);
}

void set_macro_halo_cmd(const char* macro_name,
                        int left, int bottom, int right, int top) {
  auto pine_mp = getPineMP();
  pne::Halo halo;
  halo.left = left;
  halo.bottom = bottom;
  halo.right = right;
  halo.top = top;
  pine_mp->setMacroHalo(macro_name, halo);
}

void set_soft_macros_cmd(bool enable,
                         double utilization,
                         double aspect_ratio) {
  auto pine_mp = getPineMP();
  pine_mp->enableSoftMacros(enable);
  pine_mp->setSoftMacroUtilization(utilization);
  pine_mp->setSoftMacroAspectRatio(aspect_ratio);
}

void report_soft_macros_cmd() {
  auto pine_mp = getPineMP();
  pine_mp->reportSoftMacros();
}

void set_adaptive_weighting_cmd(bool enable) {
  auto pine_mp = getPineMP();
  pine_mp->enableAdaptiveIOWeighting(enable);
}

void set_partitioning_cmd(int num_partitions,
                          double max_area_fraction,
                          int min_cells,
                          int seed,
                          bool external) {
  auto pine_mp = getPineMP();
  pine_mp->enableInternalPartitioning(!external);
  pine_mp->setPartitionTarget(num_partitions);
  pine_mp->setPartitionMaxAreaFraction(max_area_fraction);
  pine_mp->setPartitionMinCells(min_cells);
  pine_mp->setPartitionSeed(seed);
}

void report_partition_tree_cmd() {
  auto pine_mp = getPineMP();
  pine_mp->reportPartitionTree();
}

void set_corner_anchoring_cmd(bool enable) {
  auto pine_mp = getPineMP();
  pine_mp->enableCornerAnchoring(enable);
}

void set_fast_sa_cmd(bool enable, double accept_prob, double c, int k) {
  auto pine_mp = getPineMP();
  pine_mp->enableFastSA(enable);
  pine_mp->setFastSAParams(accept_prob, c, k);
}

void set_slack_moves_cmd(bool enable, double prob) {
  auto pine_mp = getPineMP();
  pine_mp->enableSlackMoves(enable);
  pine_mp->setSlackMoveProbability(prob);
}

} // namespace

%} // inline
