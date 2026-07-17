// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/RecursivePartitioner.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "par/PartitionMgr.h"
#include "utl/Logger.h"

#include "pne/SoftMacro.h"

namespace pne {

RecursivePartitioner::RecursivePartitioner(odb::dbDatabase* db,
                                           par::PartitionMgr* partition_mgr,
                                           utl::Logger* logger)
    : db_(db), partition_mgr_(partition_mgr), logger_(logger)
{
}

bool RecursivePartitioner::isSplittable(const PartitionTreeNode& node) const
{
  // Both halves must be able to hold at least min_partition_cells_.
  return node.num_cells >= std::max(2, 2 * min_partition_cells_);
}

int RecursivePartitioner::partition()
{
  tree_.clear();

  if (partition_mgr_ == nullptr) {
    logger_->warn(utl::PNE, 120, "No partitioner available for recursive bisection");
    return 0;
  }

  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    return 0;
  }
  odb::dbBlock* block = chip->getBlock();

  // Root group: every movable standard cell.
  PartitionTreeNode root;
  root.id = 0;
  for (odb::dbInst* inst : block->getInsts()) {
    if (!SoftMacroMgr::isStdCell(inst)) {
      continue;
    }
    odb::dbMaster* master = inst->getMaster();
    root.insts.push_back(inst);
    root.cell_area += static_cast<int64_t>(master->getWidth())
                      * static_cast<int64_t>(master->getHeight());
  }
  root.num_cells = static_cast<int>(root.insts.size());

  if (root.num_cells == 0) {
    logger_->warn(utl::PNE, 121, "No standard cells found for partitioning");
    return 0;
  }

  logger_->info(utl::PNE,
                122,
                "Recursive bisection: {} std cells, target partitions: {}, "
                "max partition area: {} um^2",
                root.num_cells,
                target_partitions_,
                max_partition_area_
                    / (static_cast<double>(block->getDbUnitsPerMicron())
                       * block->getDbUnitsPerMicron()));

  tree_.push_back(std::move(root));

  std::vector<int> leaves{0};
  std::unordered_set<int> unsplittable;

  // Safety cap so a mis-set size ceiling cannot recurse forever.
  constexpr int kMaxPartitions = 4096;

  while (static_cast<int>(leaves.size()) < kMaxPartitions) {
    // Pick the next leaf to split: first any oversized leaf (largest
    // first), then — while below the count target — the largest leaf.
    int candidate = -1;
    int64_t candidate_area = -1;
    bool have_oversized = false;

    for (int leaf_id : leaves) {
      const PartitionTreeNode& node = tree_[leaf_id];
      if (unsplittable.count(leaf_id) != 0 || !isSplittable(node)) {
        continue;
      }
      const bool oversized
          = max_partition_area_ > 0 && node.cell_area > max_partition_area_;
      if (have_oversized && !oversized) {
        continue;
      }
      if ((oversized && !have_oversized) || node.cell_area > candidate_area) {
        candidate = leaf_id;
        candidate_area = node.cell_area;
        have_oversized = have_oversized || oversized;
      }
    }

    const bool below_target
        = target_partitions_ > 0
          && static_cast<int>(leaves.size()) < target_partitions_;
    if (candidate < 0 || (!have_oversized && !below_target)) {
      break;
    }

    if (!bisect(candidate)) {
      unsplittable.insert(candidate);
      continue;
    }

    // Replace the split node by its children in the leaf set.
    leaves.erase(std::find(leaves.begin(), leaves.end(), candidate));
    leaves.push_back(tree_[candidate].left_child);
    leaves.push_back(tree_[candidate].right_child);
  }

  if (static_cast<int>(leaves.size()) >= kMaxPartitions) {
    logger_->warn(utl::PNE,
                  123,
                  "Recursive bisection stopped at the safety cap of {} "
                  "partitions; check the size/count limits",
                  kMaxPartitions);
  }

  // Assign final partition ids in DFS order so sibling partitions get
  // adjacent ids, then write them to the DB in the triton_part format.
  int next_partition_id = 0;
  std::vector<int> stack{0};
  while (!stack.empty()) {
    const int node_id = stack.back();
    stack.pop_back();
    PartitionTreeNode& node = tree_[node_id];
    if (node.left_child < 0 && node.right_child < 0) {
      node.leaf_partition_id = next_partition_id++;
      for (odb::dbInst* inst : node.insts) {
        if (auto* property = odb::dbIntProperty::find(inst, "partition_id")) {
          property->setValue(node.leaf_partition_id);
        } else {
          odb::dbIntProperty::create(
              inst, "partition_id", node.leaf_partition_id);
        }
      }
    } else {
      // Push right first so the left subtree is numbered first.
      if (node.right_child >= 0) {
        stack.push_back(node.right_child);
      }
      if (node.left_child >= 0) {
        stack.push_back(node.left_child);
      }
    }
  }

  int max_depth = 0;
  for (const PartitionTreeNode& node : tree_) {
    max_depth = std::max(max_depth, node.depth);
  }

  logger_->info(utl::PNE,
                124,
                "Recursive bisection created {} partitions (tree depth {})",
                next_partition_id,
                max_depth);

  return next_partition_id;
}

bool RecursivePartitioner::bisect(int node_id)
{
  PartitionTreeNode& node = tree_[node_id];
  const std::vector<odb::dbInst*>& insts = node.insts;
  const int n = static_cast<int>(insts.size());
  if (n < 2) {
    return false;
  }

  odb::dbBlock* block = db_->getChip()->getBlock();
  const double dbu = block->getDbUnitsPerMicron();

  // Local hypergraph over the group: vertices are member cells (weight =
  // area in um^2), hyperedges are nets restricted to member pins.
  // Connections leaving the group are ignored; parent-level cuts were
  // already decided at the higher level of the recursion.
  std::unordered_map<odb::dbInst*, int> inst_to_vertex;
  inst_to_vertex.reserve(n);
  std::vector<float> vertex_weights;
  vertex_weights.reserve(n);
  for (int i = 0; i < n; ++i) {
    inst_to_vertex[insts[i]] = i;
    odb::dbMaster* master = insts[i]->getMaster();
    const double area_um2 = (master->getWidth() / dbu)
                            * (master->getHeight() / dbu);
    vertex_weights.push_back(static_cast<float>(area_um2));
  }

  std::vector<std::vector<int>> hyperedges;
  std::unordered_set<odb::dbNet*> visited_nets;
  for (odb::dbInst* inst : insts) {
    for (odb::dbITerm* iterm : inst->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (net == nullptr || !visited_nets.insert(net).second) {
        continue;
      }
      const odb::dbSigType sig = net->getSigType();
      if (sig == odb::dbSigType::POWER || sig == odb::dbSigType::GROUND
          || sig == odb::dbSigType::CLOCK) {
        continue;
      }

      std::unordered_set<int> members;
      for (odb::dbITerm* net_iterm : net->getITerms()) {
        auto it = inst_to_vertex.find(net_iterm->getInst());
        if (it != inst_to_vertex.end()) {
          members.insert(it->second);
        }
        if (static_cast<int>(members.size()) >= large_net_threshold_) {
          break;
        }
      }
      if (static_cast<int>(members.size()) >= 2
          && static_cast<int>(members.size()) < large_net_threshold_) {
        hyperedges.emplace_back(members.begin(), members.end());
      }
    }
  }

  const std::vector<float> hyperedge_weights(hyperedges.size(), 1.0f);

  // Min-cut bisection with balance relaxation (same recovery scheme as
  // Hier-RTLMP): a one-sided solution means the balance constraint could
  // not be met, so loosen it and retry before falling back to a pure
  // area-balanced split.  A group without internal nets has no min-cut
  // problem to solve; split it by area directly.
  std::vector<int> solution;
  bool one_sided = true;
  float balance_constraint = 1.0f;
  constexpr float kBalanceRelaxation = 10.0f;

  while (!hyperedges.empty() && balance_constraint < 90.0f) {
    solution = partition_mgr_->PartitionKWaySimpleMode(2,
                                                       balance_constraint,
                                                       seed_,
                                                       hyperedges,
                                                       vertex_weights,
                                                       hyperedge_weights);
    int count_zero = 0;
    for (int i = 0; i < n && i < static_cast<int>(solution.size()); ++i) {
      if (solution[i] == 0) {
        ++count_zero;
      }
    }
    one_sided = static_cast<int>(solution.size()) < n || count_zero == 0
                || count_zero == n;
    if (!one_sided) {
      break;
    }
    balance_constraint += kBalanceRelaxation;
  }

  if (one_sided) {
    if (!hyperedges.empty()) {
      logger_->warn(utl::PNE,
                    125,
                    "Min-cut bisection of a {}-cell group failed to balance; "
                    "using area-balanced fallback",
                    n);
    }
    solution = fallbackSplit(insts);
  }

  // Build the two children.
  const int left_id = static_cast<int>(tree_.size());
  const int right_id = left_id + 1;

  PartitionTreeNode left;
  left.id = left_id;
  left.parent = node_id;
  left.depth = node.depth + 1;

  PartitionTreeNode right;
  right.id = right_id;
  right.parent = node_id;
  right.depth = node.depth + 1;

  for (int i = 0; i < n; ++i) {
    PartitionTreeNode& child = (solution[i] == 0) ? left : right;
    odb::dbMaster* master = insts[i]->getMaster();
    child.insts.push_back(insts[i]);
    child.cell_area += static_cast<int64_t>(master->getWidth())
                       * static_cast<int64_t>(master->getHeight());
  }
  left.num_cells = static_cast<int>(left.insts.size());
  right.num_cells = static_cast<int>(right.insts.size());

  if (left.num_cells == 0 || right.num_cells == 0) {
    return false;
  }

  // Note: `node` may dangle after push_back reallocates; use the index.
  tree_.push_back(std::move(left));
  tree_.push_back(std::move(right));
  tree_[node_id].left_child = left_id;
  tree_[node_id].right_child = right_id;
  tree_[node_id].insts.clear();
  tree_[node_id].insts.shrink_to_fit();

  return true;
}

std::vector<int> RecursivePartitioner::fallbackSplit(
    const std::vector<odb::dbInst*>& insts) const
{
  const int n = static_cast<int>(insts.size());

  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) {
    order[i] = i;
  }
  auto inst_area = [](odb::dbInst* inst) {
    odb::dbMaster* master = inst->getMaster();
    return static_cast<int64_t>(master->getWidth())
           * static_cast<int64_t>(master->getHeight());
  };
  std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
    return inst_area(insts[lhs]) > inst_area(insts[rhs]);
  });

  // Greedy: biggest cells first, always into the lighter half.
  std::vector<int> solution(n, 0);
  int64_t area[2] = {0, 0};
  for (int idx : order) {
    const int side = (area[0] <= area[1]) ? 0 : 1;
    solution[idx] = side;
    area[side] += inst_area(insts[idx]);
  }
  return solution;
}

void reportPartitionTree(const std::vector<PartitionTreeNode>& tree,
                         utl::Logger* logger)
{
  if (tree.empty()) {
    logger->warn(utl::PNE, 126,
                 "No partition tree available. "
                 "Run pine_mp with soft macros enabled first.");
    return;
  }

  logger->info(utl::PNE, 127, "Partition tree ({} nodes):", tree.size());

  // DFS from the root, printing depth-indented nodes.
  std::vector<int> stack{0};
  while (!stack.empty()) {
    const int node_id = stack.back();
    stack.pop_back();
    const PartitionTreeNode& node = tree[node_id];

    const std::string indent(2 * node.depth, ' ');
    if (node.left_child < 0 && node.right_child < 0) {
      logger->info(utl::PNE,
                   128,
                   "{}[leaf] partition {}: {} cells, area {} DBU^2",
                   indent,
                   node.leaf_partition_id,
                   node.num_cells,
                   node.cell_area);
    } else {
      logger->info(utl::PNE,
                   129,
                   "{}[node {}] depth {}: {} cells, area {} DBU^2",
                   indent,
                   node.id,
                   node.depth,
                   node.num_cells,
                   node.cell_area);
      if (node.right_child >= 0) {
        stack.push_back(node.right_child);
      }
      if (node.left_child >= 0) {
        stack.push_back(node.left_child);
      }
    }
  }
}

}  // namespace pne
