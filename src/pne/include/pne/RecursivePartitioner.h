// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <vector>

#include "odb/db.h"

namespace par {
class PartitionMgr;
}

namespace utl {
class Logger;
}

namespace pne {

// One node of the recursive-bisection hierarchy.  Leaves are the final
// partitions (consumed as soft macros); internal nodes record how the
// netlist was divided, preserved for future hierarchical placement.
struct PartitionTreeNode
{
  int id = -1;
  int parent = -1;
  int depth = 0;
  int left_child = -1;   // -1 on leaves
  int right_child = -1;
  int leaf_partition_id = -1;  // final partition id; -1 on internal nodes
  int num_cells = 0;
  int64_t cell_area = 0;  // sum of member master areas (DBU^2)
  // Member instances.  Kept on leaves; cleared once a node is split.
  std::vector<odb::dbInst*> insts;
};

// Log the bisection hierarchy (indented by depth).
void reportPartitionTree(const std::vector<PartitionTreeNode>& tree,
                         utl::Logger* logger);

// Recursive bisection of the standard-cell netlist using the par module's
// in-memory 2-way partitioner (the same engine behind triton_part).  This
// replaces the external triton_part call in the flow: starting from all
// std cells, the largest group is repeatedly split in half until the
// requested number of partitions and/or the maximum partition size is
// reached.  Divide and conquer is faster than one flat k-way call and the
// bisection tree is retained as a hierarchy for future hierarchical
// placement.
//
// The final partition ids are written as "partition_id" dbIntProperty on
// every member instance — the same format triton_part produces and
// SoftMacroMgr::buildFromPartitions consumes.
class RecursivePartitioner
{
 public:
  RecursivePartitioner(odb::dbDatabase* db,
                       par::PartitionMgr* partition_mgr,
                       utl::Logger* logger);

  // Split until at least `num` leaves exist (0 disables the count target).
  void setTargetPartitions(int num) { target_partitions_ = num; }
  // Keep splitting any partition whose total cell area exceeds this value
  // (DBU^2; 0 disables the ceiling).  Takes precedence over the count
  // target: oversized partitions are split even beyond it.
  void setMaxPartitionArea(int64_t area) { max_partition_area_ = area; }
  // Never split groups below this cell count.
  void setMinPartitionCells(int num) { min_partition_cells_ = num; }
  void setSeed(unsigned seed) { seed_ = seed; }
  // Nets with this many or more member pins inside a group are ignored
  // when building the bisection hypergraph (clock/reset-like fanout
  // carries no locality information).
  void setLargeNetThreshold(int threshold) { large_net_threshold_ = threshold; }

  // Run the recursive bisection over all movable std cells.
  // Returns the number of partitions created (0 on failure).
  int partition();

  const std::vector<PartitionTreeNode>& getTree() const { return tree_; }

 private:
  bool isSplittable(const PartitionTreeNode& node) const;
  // Split tree_[node_id] in two; returns false when no usable split exists.
  bool bisect(int node_id);
  // Area-balanced greedy halving used when the min-cut partitioner returns
  // a one-sided solution even at the loosest balance constraint.
  std::vector<int> fallbackSplit(const std::vector<odb::dbInst*>& insts) const;

  odb::dbDatabase* db_;
  par::PartitionMgr* partition_mgr_;
  utl::Logger* logger_;

  int target_partitions_ = 0;
  int64_t max_partition_area_ = 0;
  int min_partition_cells_ = 50;
  unsigned seed_ = 1;
  int large_net_threshold_ = 100;

  std::vector<PartitionTreeNode> tree_;
};

}  // namespace pne
