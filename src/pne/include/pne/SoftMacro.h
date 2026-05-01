// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "odb/db.h"
#include "utl/Logger.h"

namespace pne {

// SoftMacro represents a group of standard cells that belong to the same
// par partition. It acts as a virtual (synthetic) macro in the B*-Tree
// placement engine so that clusters of stdcells can be co-placed with hard
// macros during floorplanning.
//
// Dimensions are derived from the cumulative area of member cells and a
// user-controlled target utilization and aspect ratio.  The actual physical
// placement of individual cells within the soft macro region is left to
// later global / detailed placement steps.
struct SoftMacro
{
  int partition_id = -1;

  // Standard-cell instances that belong to this partition.
  // Hard macros are excluded; they are handled by BStarTree as usual.
  std::vector<odb::dbInst*> instances;

  // Soft macro bounding box dimensions (in DBU) used during placement.
  // Set by SoftMacroMgr::buildFromPartitions().
  int width = 0;
  int height = 0;

  // Placed lower-left origin (in DBU).  Written by BStarTree::applyPlacement()
  // so callers can query the final position of each soft macro.
  int x = 0;
  int y = 0;

  // Human-readable identifier.
  std::string getName() const
  {
    return "soft_macro_" + std::to_string(partition_id);
  }

  // Sum of the master bounding-box areas of all member instances (DBU^2).
  int64_t getCellArea() const;

  // Placed area (width * height) after dimensions have been computed.
  int64_t getArea() const
  {
    return static_cast<int64_t>(width) * static_cast<int64_t>(height);
  }
};

// SoftMacroMgr constructs SoftMacro objects from the "partition_id"
// dbIntProperty that triton_part (par module) writes on every dbInst after a
// successful partitioning run.  Only standard-cell instances (non-macro,
// non-pad, non-cover) are considered; hard macros keep their own tree nodes.
//
// Typical usage:
//   SoftMacroMgr mgr(db, logger);
//   int n = mgr.buildFromPartitions(0.7, 1.0);
//   if (n > 0)
//     // attach soft macros to the B*-Tree ...
class SoftMacroMgr
{
 public:
  SoftMacroMgr(odb::dbDatabase* db, utl::Logger* logger);

  // Build soft macros from existing "partition_id" properties in the DB.
  //
  // target_utilization  – target stdcell packing density inside the soft macro
  //                        rectangle (0 < u <= 1, default 0.7).
  // aspect_ratio        – width / height of the soft macro rectangle
  //                        (default 1.0 → square).
  //
  // Returns the number of soft macros created; 0 if no "partition_id"
  // properties are found (i.e. triton_part was not run first).
  int buildFromPartitions(double target_utilization, double aspect_ratio);

  const std::vector<SoftMacro>& getSoftMacros() const { return soft_macros_; }
  std::vector<SoftMacro>& getSoftMacros() { return soft_macros_; }

  bool hasSoftMacros() const { return !soft_macros_.empty(); }
  int getNumSoftMacros() const
  {
    return static_cast<int>(soft_macros_.size());
  }

  // Print per-partition statistics via the logger (INFO level).
  void reportStats() const;

 private:
  odb::dbDatabase* db_;
  utl::Logger* logger_;
  std::vector<SoftMacro> soft_macros_;

  // Returns true for instances that should be grouped into soft macros.
  // Hard macros, pads, covers, and rings are excluded.
  static bool isStdCell(odb::dbInst* inst);

  // Compute (width, height) from cumulative cell area, target utilization,
  // and desired aspect ratio.
  static std::pair<int, int> computeDimensions(int64_t cell_area,
                                               double utilization,
                                               double aspect_ratio);
};

}  // namespace pne
