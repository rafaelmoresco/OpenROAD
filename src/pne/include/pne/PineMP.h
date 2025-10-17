// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include "odb/db.h"
#include "utl/Logger.h"

namespace pne {

class PineMP
{
public:
  PineMP(utl::Logger* logger);
  ~PineMP() = default;

  // Basic placement function - to be implemented
  bool place(int num_threads = 1);

  // Utility functions
  void setDebug(bool debug) { debug_ = debug; }
  bool getDebug() const { return debug_; }

private:
  utl::Logger* logger_;
  bool debug_ = false;
};

}  // namespace pne
