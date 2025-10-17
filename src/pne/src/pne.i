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

} // namespace

%} // inline
