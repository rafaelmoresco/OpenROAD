// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include "odb/db.h"
#include "utl/Logger.h"

namespace odb {
class dbDatabase;
class dbInst;
class dbOrientType;
}  // namespace odb

namespace sta {
class dbNetwork;
class dbSta;
}  // namespace sta

namespace par {
class PartitionMgr;
}
namespace pne {

class PineMP
{
public:
  PineMP(odb::dbDatabase* db, utl::Logger* logger);
  ~PineMP() = default;

  // Basic placement function - to be implemented
  bool place(int num_threads = 1);

  // Utility functions
  void setDebug(bool debug) { debug_ = debug; }
  bool getDebug() const { return debug_; }

private:
  utl::Logger* logger_;
  odb::dbDatabase* db_;
  bool debug_ = false;
};

}  // namespace pne
