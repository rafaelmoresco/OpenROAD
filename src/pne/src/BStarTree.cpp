// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/BStarTree.h"

#include <algorithm>
#include <cstdint>
#include <limits>

#include "pne/SoftMacro.h"

namespace pne {

namespace {

bool isDescendant(BStarNode* root, BStarNode* target)
{
  if (root == nullptr || target == nullptr) {
    return false;
  }
  if (root == target) {
    return true;
  }
  return isDescendant(root->getLeft(), target)
         || isDescendant(root->getRight(), target);
}

BStarNode* rightmostNode(BStarNode* node)
{
  if (node == nullptr) {
    return nullptr;
  }
  while (node->getRight() != nullptr) {
    node = node->getRight();
  }
  return node;
}

void attachAsRightmost(BStarNode* host, BStarNode* child)
{
  if (host == nullptr || child == nullptr) {
    return;
  }

  BStarNode* tail = rightmostNode(host);
  tail->setRight(child);
  child->setParent(tail);
}

}  // namespace

//-----------------------------------------------------------------------------
// BStarNode implementation
//-----------------------------------------------------------------------------

BStarNode::BStarNode(odb::dbInst* inst, int id)
    : id_(id), inst_(inst), local_orient_(odb::dbOrientType::R0)
{
}

BStarNode::BStarNode(SoftMacro* sm, int id)
    : id_(id), inst_(nullptr), soft_macro_(sm), local_orient_(odb::dbOrientType::R0)
{
}

int BStarNode::getMacroWidth() const
{
  if (soft_macro_ != nullptr) {
    return soft_macro_->width;
  }
  odb::dbBox* bbox = inst_->getBBox();
  return bbox->getBox().dx();
}

int BStarNode::getMacroHeight() const
{
  if (soft_macro_ != nullptr) {
    return soft_macro_->height;
  }
  odb::dbBox* bbox = inst_->getBBox();
  return bbox->getBox().dy();
}

int BStarNode::getWidth() const
{
  return getMacroWidth() + halo_.left + halo_.right;
}

int BStarNode::getHeight() const
{
  return getMacroHeight() + halo_.bottom + halo_.top;
}

odb::dbOrientType BStarNode::getOrientation() const
{
  if (soft_macro_ != nullptr) {
    return local_orient_;
  }
  return inst_->getOrient();
}

void BStarNode::setOrientation(odb::dbOrientType orient)
{
  if (soft_macro_ != nullptr) {
    // Soft macros have no physical orientation in the DB.
    // Store locally so save/restore round-trips are consistent.
    local_orient_ = orient;
    return;
  }
  inst_->setOrient(orient);
}

//-----------------------------------------------------------------------------
// BStarTree implementation
//-----------------------------------------------------------------------------

BStarTree::BStarTree()
{
}

BStarTree::~BStarTree()
{
  clear();
}

void BStarTree::clear()
{
  nodes_.clear();
  root_ = nullptr;
  width_ = 0;
  height_ = 0;
  contours_.clear();
}

void BStarTree::addMacro(odb::dbInst* inst)
{
  int id = nodes_.size();
  auto node = std::make_unique<BStarNode>(inst, id);

  if (root_ == nullptr) {
    root_ = node.get();
  } else {
    // Heap-based balanced binary tree: parent of node[i] is node[(i-1)/2].
    // Node[i] is left child  when i is odd, right child when i is even.
    int parent_id = (id - 1) / 2;
    BStarNode* parent = findNode(parent_id);
    if (id % 2 == 1) {
      parent->setLeft(node.get());
    } else {
      parent->setRight(node.get());
    }
    node->setParent(parent);
  }

  applyDefaultHalo(node.get());
  nodes_.push_back(std::move(node));
}

void BStarTree::addSoftMacro(SoftMacro* sm)
{
  int id = nodes_.size();
  auto node = std::make_unique<BStarNode>(sm, id);

  if (root_ == nullptr) {
    root_ = node.get();
  } else {
    int parent_id = (id - 1) / 2;
    BStarNode* parent = findNode(parent_id);
    if (id % 2 == 1) {
      parent->setLeft(node.get());
    } else {
      parent->setRight(node.get());
    }
    node->setParent(parent);
  }

  applyDefaultHalo(node.get());
  nodes_.push_back(std::move(node));
}

void BStarTree::addSoftMacros(const std::vector<SoftMacro*>& soft_macros,
                              int die_height,
                              int halo_y)
{
  if (soft_macros.empty()) {
    return;
  }

  BStarNode* attach_point = nullptr;
  if (root_ != nullptr) {
    // Find the rightmost node in the existing left-child chain. Soft macros
    // should extend the existing layout to the right in a structured way.
    attach_point = root_;
    while (attach_point->getLeft() != nullptr) {
      attach_point = attach_point->getLeft();
    }
  }

  int min_h = std::numeric_limits<int>::max();
  for (SoftMacro* sm : soft_macros) {
    min_h = std::min(min_h, sm->height);
  }

  const int tall_threshold = die_height - min_h;
  std::vector<SoftMacro*> tall_macros;
  std::vector<SoftMacro*> stack_macros;
  for (SoftMacro* sm : soft_macros) {
    if (sm->height > tall_threshold) {
      tall_macros.push_back(sm);
    } else {
      stack_macros.push_back(sm);
    }
  }

  for (SoftMacro* sm : tall_macros) {
    int id = static_cast<int>(nodes_.size());
    auto node = std::make_unique<BStarNode>(sm, id);
    if (attach_point == nullptr) {
      root_ = node.get();
    } else {
      attach_point->setLeft(node.get());
      node->setParent(attach_point);
    }
    attach_point = node.get();
    nodes_.push_back(std::move(node));
  }

  if (stack_macros.empty()) {
    return;
  }

  int max_stack_h = 0;
  for (SoftMacro* sm : stack_macros) {
    max_stack_h = std::max(max_stack_h, sm->height);
  }
  const int effective_slot_h = max_stack_h + 2 * halo_y;
  const int max_depth = std::max(1, die_height / effective_slot_h);
  const int n = static_cast<int>(stack_macros.size());
  const int n_cols = (n + max_depth - 1) / max_depth;

  BStarNode* prev_col_head = nullptr;
  for (int col = 0; col < n_cols; col++) {
    const int start = col * max_depth;
    const int end = std::min(n, start + max_depth);

    int id = static_cast<int>(nodes_.size());
    auto head_node = std::make_unique<BStarNode>(stack_macros[start], id);
    BStarNode* col_head = head_node.get();

    if (prev_col_head == nullptr) {
      if (attach_point == nullptr) {
        root_ = col_head;
      } else {
        attach_point->setLeft(col_head);
        col_head->setParent(attach_point);
      }
    } else {
      prev_col_head->setLeft(col_head);
      col_head->setParent(prev_col_head);
    }

    nodes_.push_back(std::move(head_node));
    applyDefaultHalo(nodes_.back().get());
    prev_col_head = col_head;

    BStarNode* cur = col_head;
    for (int i = start + 1; i < end; i++) {
      int nid = static_cast<int>(nodes_.size());
      auto node = std::make_unique<BStarNode>(stack_macros[i], nid);
      cur->setRight(node.get());
      node->setParent(cur);
      cur = node.get();
      nodes_.push_back(std::move(node));
      applyDefaultHalo(nodes_.back().get());
    }
  }
}

void BStarTree::buildFromMacros(const std::vector<odb::dbInst*>& macros,
                                 int die_height,
                                 int halo_y)
{
  clear();
  if (macros.empty()) {
    return;
  }

  // Compute minimum macro height to determine what can be stacked.
  int min_h = std::numeric_limits<int>::max();
  for (auto* inst : macros) {
    min_h = std::min(min_h, inst->getBBox()->getBox().dy());
  }

  // A "tall" macro is one where even the shortest sibling macro cannot be
  // stacked above it within the die (no room left after placing the tall macro).
  const int tall_threshold = die_height - min_h;

  std::vector<odb::dbInst*> tall_macros;
  std::vector<odb::dbInst*> stack_macros;
  for (auto* inst : macros) {
    if (inst->getBBox()->getBox().dy() > tall_threshold) {
      tall_macros.push_back(inst);
    } else {
      stack_macros.push_back(inst);
    }
  }

  // Phase 1: tall macros as a pure left-child horizontal chain.
  // Left children in B*-tree go to larger X; right children go above (same X).
  // Tall macros MUST NOT have right-child descendants (would overflow die height).
  for (auto* inst : tall_macros) {
    int id = static_cast<int>(nodes_.size());
    auto node = std::make_unique<BStarNode>(inst, id);
    if (root_ == nullptr) {
      root_ = node.get();
    } else {
      BStarNode* tail = root_;
      while (tail->getLeft() != nullptr) {
        tail = tail->getLeft();
      }
      tail->setLeft(node.get());
      node->setParent(tail);
    }
    applyDefaultHalo(node.get());
    nodes_.push_back(std::move(node));
  }

  if (stack_macros.empty()) {
    return;
  }

  // Phase 2: stackable macros in vertical columns (right-child chains).
  // Compute how many macros fit in one column, accounting for the per-side
  // vertical halo (top + bottom = 2*halo_y) that will be added to each macro
  // after the tree is built.  Without this correction the initial packed height
  // can already exceed the core boundary before SA starts.
  int max_stack_h = 0;
  for (auto* inst : stack_macros) {
    max_stack_h = std::max(max_stack_h, inst->getBBox()->getBox().dy());
  }
  const int effective_slot_h = max_stack_h + 2 * halo_y;
  const int max_depth = std::max(1, die_height / effective_slot_h);
  const int n = static_cast<int>(stack_macros.size());
  const int n_cols = (n + max_depth - 1) / max_depth;

  // Each column head is a left child of the previous column head (or of the
  // last tall macro if this is the first column).  Within a column the macros
  // form a right-child chain starting from the column head.
  BStarNode* prev_col_head = nullptr;

  for (int col = 0; col < n_cols; col++) {
    const int start = col * max_depth;
    const int end = std::min(n, start + max_depth);

    // Create column head.
    int id = static_cast<int>(nodes_.size());
    auto head_node = std::make_unique<BStarNode>(stack_macros[start], id);
    BStarNode* col_head = head_node.get();

    if (prev_col_head == nullptr) {
      // First column: attach as left child of the last node in the tall chain,
      // or as root if there were no tall macros.
      if (root_ == nullptr) {
        root_ = col_head;
      } else {
        BStarNode* tail = root_;
        while (tail->getLeft() != nullptr) {
          tail = tail->getLeft();
        }
        tail->setLeft(col_head);
        col_head->setParent(tail);
      }
    } else {
      prev_col_head->setLeft(col_head);
      col_head->setParent(prev_col_head);
    }

    nodes_.push_back(std::move(head_node));
    applyDefaultHalo(nodes_.back().get());
    prev_col_head = col_head;

    // Remaining macros in this column as a right-child chain.
    BStarNode* cur = col_head;
    for (int i = start + 1; i < end; i++) {
      int nid = static_cast<int>(nodes_.size());
      auto node = std::make_unique<BStarNode>(stack_macros[i], nid);
      cur->setRight(node.get());
      node->setParent(cur);
      cur = node.get();
      nodes_.push_back(std::move(node));
      applyDefaultHalo(nodes_.back().get());
    }
  }
}

void BStarTree::pack()
{
  if (root_ == nullptr) {
    width_ = 0;
    height_ = 0;
    return;
  }
  
  // Reset contours
  contours_.clear();
  contours_.push_back({0, 0, std::numeric_limits<int>::max()});

  std::vector<bool> visited(nodes_.size(), false);
  
  // Pack starting from root at origin
  packRecursive(root_, 0, visited);

  // If the tree was corrupted by a perturbation, pack any disconnected nodes
  // as separate roots to avoid leaving stale coordinates.
  for (const auto& node : nodes_) {
    const int id = node->getId();
    if (id >= 0 && id < static_cast<int>(visited.size()) && !visited[id]) {
      packRecursive(node.get(), 0, visited);
    }
  }
  
  // Compute bounding box
  width_ = 0;
  height_ = 0;
  for (const auto& node : nodes_) {
    const int64_t right64 = static_cast<int64_t>(node->getX())
                            + static_cast<int64_t>(node->getWidth());
    const int64_t top64 = static_cast<int64_t>(node->getY())
                          + static_cast<int64_t>(node->getHeight());
    const int right = static_cast<int>(std::max<int64_t>(
        0,
        std::min<int64_t>(right64, std::numeric_limits<int>::max())));
    const int top = static_cast<int>(std::max<int64_t>(
        0,
        std::min<int64_t>(top64, std::numeric_limits<int>::max())));
    width_ = std::max(width_, right);
    height_ = std::max(height_, top);
  }
}

void BStarTree::packRecursive(BStarNode* node, int x, std::vector<bool>& visited)
{
  if (node == nullptr) {
    return;
  }

  const int id = node->getId();
  if (id < 0 || id >= static_cast<int>(visited.size())) {
    return;
  }
  if (visited[id]) {
    return;
  }
  visited[id] = true;
  
  // Find Y coordinate using contour
  int y = findContourY(x, node->getWidth());
  node->setX(x);
  node->setY(y);
  
  // Update contour
  updateContour(x, y, node->getWidth(), node->getHeight());
  
  // Left child: place to the right (same y-level)
  if (node->getLeft() != nullptr) {
    const int64_t child_x64 = static_cast<int64_t>(x)
                              + static_cast<int64_t>(node->getWidth());
    const int child_x = static_cast<int>(std::max<int64_t>(
        std::numeric_limits<int>::min(),
        std::min<int64_t>(child_x64, std::numeric_limits<int>::max())));
    packRecursive(node->getLeft(), child_x, visited);
  }
  
  // Right child: place above (compacted upward)
  if (node->getRight() != nullptr) {
    packRecursive(node->getRight(), x, visited);
  }
}

int BStarTree::findContourY(int x, int width)
{
  int max_y = 0;
  const int64_t range_begin = static_cast<int64_t>(x);
  const int64_t range_end = static_cast<int64_t>(x) + static_cast<int64_t>(width);
  
  // Find maximum Y in the range [x, x+width)
  for (const auto& c : contours_) {
    if (c.width <= 0) {
      continue;
    }
    const int64_t c_begin = static_cast<int64_t>(c.x);
    const int64_t c_end = static_cast<int64_t>(c.x)
                          + static_cast<int64_t>(c.width);
    
    // Check if this contour segment overlaps with [x, x+width)
    if (c_begin < range_end && c_end > range_begin) {
      max_y = std::max(max_y, c.y);
    }
  }
  
  return max_y;
}

void BStarTree::updateContour(int x, int y, int width, int height)
{
  std::vector<Contour> new_contours;
  int new_y = y + height;
  const int64_t block_begin = static_cast<int64_t>(x);
  const int64_t block_end = static_cast<int64_t>(x) + static_cast<int64_t>(width);
  
  for (const auto& c : contours_) {
    if (c.width <= 0) {
      continue;
    }
    const int64_t c_begin = static_cast<int64_t>(c.x);
    const int64_t c_end = static_cast<int64_t>(c.x)
                          + static_cast<int64_t>(c.width);
    
    // Contour completely before block
    if (c_end <= block_begin) {
      new_contours.push_back(c);
    }
    // Contour completely after block
    else if (c_begin >= block_end) {
      new_contours.push_back(c);
    }
    // Contour overlaps with block
    else {
      // Keep part before block
      if (c_begin < block_begin) {
        const int64_t before_width = block_begin - c_begin;
        if (before_width > 0) {
          new_contours.push_back(
              {c.x, c.y, static_cast<int>(std::min<int64_t>(
                            before_width,
                            std::numeric_limits<int>::max()))});
        }
      }
      
      // Keep part after block
      if (c_end > block_end) {
        const int64_t after_width = c_end - block_end;
        if (after_width > 0) {
          const int block_start = static_cast<int>(std::max<int64_t>(
              std::numeric_limits<int>::min(),
              std::min<int64_t>(block_end, std::numeric_limits<int>::max())));
          new_contours.push_back(
              {block_start, c.y, static_cast<int>(std::min<int64_t>(
                                   after_width,
                                   std::numeric_limits<int>::max()))});
        }
      }
    }
  }
  
  // Add new contour for the block
  new_contours.push_back({x, new_y, width});
  
  // Sort by x coordinate
  std::sort(new_contours.begin(), new_contours.end(),
            [](const Contour& a, const Contour& b) { return a.x < b.x; });
  
  // Merge adjacent contours with same height
  contours_.clear();
  for (const auto& c : new_contours) {
    if (!contours_.empty() && 
        contours_.back().y == c.y && 
        contours_.back().x + contours_.back().width == c.x) {
      contours_.back().width += c.width;
    } else {
      contours_.push_back(c);
    }
  }
}

BStarNode* BStarTree::findNode(int id)
{
  if (id >= 0 && id < static_cast<int>(nodes_.size())) {
    return nodes_[id].get();
  }
  return nullptr;
}

void BStarTree::swapNodes(int id1, int id2)
{
  BStarNode* node1 = findNode(id1);
  BStarNode* node2 = findNode(id2);

  if (node1 == nullptr || node2 == nullptr || node1 == node2) {
    return;
  }

  // Walk parent chains to detect ancestor-descendant relationship.
  // Swapping a node with its own ancestor/descendant would create a cycle.
  auto isAncOf = [](BStarNode* anc, BStarNode* node) {
    BStarNode* cur = node->getParent();
    int count = 0;
    while (cur) {
      if (++count > 100000) {
        // In case the tree is already broken
        return true;
      }
      if (cur == anc) {
        return true;
      }
      cur = cur->getParent();
    }
    return false;
  };
  if (isAncOf(node1, node2) || isAncOf(node2, node1)) {
    return;
  }

  BStarNode* parent1 = node1->getParent();
  BStarNode* parent2 = node2->getParent();
  const bool isLeft1 = parent1 && (parent1->getLeft() == node1);
  const bool isLeft2 = parent2 && (parent2->getLeft() == node2);

  // Capture children before we start rewiring.
  BStarNode* left1  = node1->getLeft();
  BStarNode* right1 = node1->getRight();
  BStarNode* left2  = node2->getLeft();
  BStarNode* right2 = node2->getRight();

  // Place node2 in node1's old slot.
  node2->setParent(parent1);
  if (parent1) {
    if (isLeft1) {
      parent1->setLeft(node2);
    } else {
      parent1->setRight(node2);
    }
  } else {
    root_ = node2;
  }

  // Place node1 in node2's old slot.
  node1->setParent(parent2);
  if (parent2) {
    if (isLeft2) {
      parent2->setLeft(node1);
    } else {
      parent2->setRight(node1);
    }
  } else {
    root_ = node1;
  }

  // Exchange the two subtrees by swapping children.
  node1->setLeft(left2);
  node1->setRight(right2);
  node2->setLeft(left1);
  node2->setRight(right1);

  // Fix children's parent pointers.
  if (node1->getLeft())  { node1->getLeft()->setParent(node1);  }
  if (node1->getRight()) { node1->getRight()->setParent(node1); }
  if (node2->getLeft())  { node2->getLeft()->setParent(node2);  }
  if (node2->getRight()) { node2->getRight()->setParent(node2); }
}

void BStarTree::rotateNode(int id)
{
  BStarNode* node = findNode(id);
  if (node == nullptr) {
    return;
  }
  
  // Simple rotation: flip orientation 90 degrees
  odb::dbOrientType current = node->getOrientation();
  odb::dbOrientType new_orient;
  
  switch (current.getValue()) {
    case odb::dbOrientType::R0:
      new_orient = odb::dbOrientType::R90;
      break;
    case odb::dbOrientType::R90:
      new_orient = odb::dbOrientType::R180;
      break;
    case odb::dbOrientType::R180:
      new_orient = odb::dbOrientType::R270;
      break;
    case odb::dbOrientType::R270:
      new_orient = odb::dbOrientType::R0;
      break;
    default:
      new_orient = odb::dbOrientType::R0;
  }
  
  node->setOrientation(new_orient);
}

void BStarTree::removeFromTree(BStarNode* node)
{
  if (node == nullptr) {
    return;
  }
  
  BStarNode* parent = node->getParent();
  BStarNode* left = node->getLeft();
  BStarNode* right = node->getRight();
  
  // Update parent's pointer
  if (parent) {
    if (parent->getLeft() == node) {
      parent->setLeft(nullptr);
    } else {
      parent->setRight(nullptr);
    }
  }
  
  node->setParent(nullptr);
  node->setLeft(nullptr);
  node->setRight(nullptr);
  
  // Reattach children (simplified strategy)
  if (parent) {
    if (left) {
      if (parent->getLeft() == nullptr) {
        parent->setLeft(left);
      } else {
        parent->setRight(left);
      }
      left->setParent(parent);
    }
    if (right) {
      if (parent->getRight() == nullptr) {
        parent->setRight(right);
      } else {
        // Find a place in tree
        BStarNode* current = root_;
        while (current->getRight() != nullptr) {
          current = current->getRight();
        }
        current->setRight(right);
      }
      right->setParent(parent);
    }
  }
}

void BStarTree::insertInTree(BStarNode* node, BStarNode* parent, bool as_left)
{
  if (node == nullptr) {
    return;
  }
  
  if (parent == nullptr) {
    root_ = node;
    node->setParent(nullptr);
    return;
  }
  
  if (as_left) {
    BStarNode* old_left = parent->getLeft();
    parent->setLeft(node);
    node->setParent(parent);
    if (old_left) {
      node->setLeft(old_left);
      old_left->setParent(node);
    }
  } else {
    BStarNode* old_right = parent->getRight();
    parent->setRight(node);
    node->setParent(parent);
    if (old_right) {
      node->setRight(old_right);
      old_right->setParent(node);
    }
  }
}

void BStarTree::moveNode(int id, int new_parent_id, bool as_left_child)
{
  BStarNode* node = findNode(id);
  BStarNode* new_parent = findNode(new_parent_id);
  
  if (node == nullptr || new_parent == nullptr) {
    return;
  }
  
  // Cannot move to self or under own subtree.
  if (node == new_parent) {
    return;
  }

  if (isDescendant(node, new_parent)) {
    return;
  }

  BStarNode* old_parent = node->getParent();
  BStarNode* old_left = node->getLeft();
  BStarNode* old_right = node->getRight();

  // Detach node and stitch its former location with one child to
  // keep the remaining tree connected.
  BStarNode* replacement = old_left;
  if (replacement != nullptr) {
    replacement->setParent(old_parent);
    attachAsRightmost(replacement, old_right);
  } else {
    replacement = old_right;
    if (replacement != nullptr) {
      replacement->setParent(old_parent);
    }
  }

  if (old_parent != nullptr) {
    if (old_parent->getLeft() == node) {
      old_parent->setLeft(replacement);
    } else if (old_parent->getRight() == node) {
      old_parent->setRight(replacement);
    }
  } else {
    root_ = replacement;
  }

  node->setParent(nullptr);
  node->setLeft(nullptr);
  node->setRight(nullptr);

  // Insert under new parent and preserve any displaced subtree.
  BStarNode* displaced = as_left_child ? new_parent->getLeft()
                                       : new_parent->getRight();
  if (as_left_child) {
    new_parent->setLeft(node);
  } else {
    new_parent->setRight(node);
  }
  node->setParent(new_parent);

  if (displaced != nullptr && displaced != node) {
    node->setRight(displaced);
    displaced->setParent(node);
  }
  
  if (root_ == nullptr) {
    root_ = node;
    node->setParent(nullptr);
  }
}

void BStarTree::applyPlacement()
{
  for (const auto& node : nodes_) {
    // The packed coordinates include halo padding.
    // The actual macro origin is offset by the left/bottom halo.
    const Halo& halo = node->getHalo();
    const int x = node->getX() + halo.left;
    const int y = node->getY() + halo.bottom;

    if (node->isSoftMacro()) {
      // Store the placed origin in the SoftMacro so callers can query it.
      // Individual stdcells are placed within this region by global placement.
      node->getSoftMacro()->x = x;
      node->getSoftMacro()->y = y;
    } else {
      node->getInst()->setLocation(x, y);
      node->getInst()->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    }
  }
}

void BStarTree::save()
{
  saveSnapshot(SnapshotSlot::CURRENT);
}

void BStarTree::restore()
{
  restoreSnapshot(SnapshotSlot::CURRENT);
}

void BStarTree::saveSnapshot(SnapshotSlot slot)
{
  auto& backup = backups_[static_cast<size_t>(slot)];

  backup.clear();
  backup.reserve(nodes_.size());

  for (const auto& node : nodes_) {
    NodeState state;
    state.left_id = node->getLeft() ? node->getLeft()->getId() : -1;
    state.right_id = node->getRight() ? node->getRight()->getId() : -1;
    state.parent_id = node->getParent() ? node->getParent()->getId() : -1;
    state.x = node->getX();
    state.y = node->getY();
    state.orient = node->getOrientation();
    state.halo = node->getHalo();
    backup.push_back(state);
  }
}

void BStarTree::restoreSnapshot(SnapshotSlot slot)
{
  const auto& backup = backups_[static_cast<size_t>(slot)];

  if (backup.size() != nodes_.size()) {
    return;
  }
  
  // First pass: restore orientations, coordinates, and halos
  for (size_t i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->setX(backup[i].x);
    nodes_[i]->setY(backup[i].y);
    nodes_[i]->setOrientation(backup[i].orient);
    nodes_[i]->setHalo(backup[i].halo);
  }
  
  // Second pass: restore tree structure
  for (size_t i = 0; i < nodes_.size(); ++i) {
    BStarNode* left = findNode(backup[i].left_id);
    BStarNode* right = findNode(backup[i].right_id);
    BStarNode* parent = findNode(backup[i].parent_id);
    
    nodes_[i]->setLeft(left);
    nodes_[i]->setRight(right);
    nodes_[i]->setParent(parent);
    
    if (backup[i].parent_id == -1) {
      root_ = nodes_[i].get();
    }
  }
}

// Determine which sides of a macro instance have signal pins.
// Returns a Halo where each direction is either halo_value (if pins
// exist on that side) or 0.
static Halo computeSideHalo(odb::dbInst* inst, int halo_x, int halo_y)
{
  odb::Rect master_box;
  inst->getMaster()->getPlacementBoundary(master_box);
  if (master_box.dx() == 0 || master_box.dy() == 0) {
    odb::dbMaster* master = inst->getMaster();
    master_box = odb::Rect(0, 0, master->getWidth(), master->getHeight());
  }

  const int mw = master_box.dx();
  const int mh = master_box.dy();

  bool has_left = false;
  bool has_right = false;
  bool has_bottom = false;
  bool has_top = false;

  for (odb::dbITerm* iterm : inst->getITerms()) {
    // Only consider signal pins for halo
    const odb::dbSigType sig = iterm->getSigType();
    if (sig == odb::dbSigType::POWER || sig == odb::dbSigType::GROUND) {
      continue;
    }

    odb::dbMTerm* mterm = iterm->getMTerm();
    if (mterm == nullptr) {
      continue;
    }

    // Compute pin center in master coordinates
    odb::Rect pin_bbox = mterm->getBBox();
    int cx = (pin_bbox.xMin() + pin_bbox.xMax()) / 2;
    int cy = (pin_bbox.yMin() + pin_bbox.yMax()) / 2;

    // Map to fractional position [0..1] within the master cell
    double fx = (mw > 0) ? static_cast<double>(cx - master_box.xMin()) / mw
                         : 0.5;
    double fy = (mh > 0) ? static_cast<double>(cy - master_box.yMin()) / mh
                         : 0.5;

    // Determine dominant side(s) — a pin near a boundary edge contributes
    // to that side.  Use a threshold of 15% of the cell dimension.
    constexpr double kEdgeThreshold = 0.15;

    if (fx < kEdgeThreshold) {
      has_left = true;
    }
    if (fx > (1.0 - kEdgeThreshold)) {
      has_right = true;
    }
    if (fy < kEdgeThreshold) {
      has_bottom = true;
    }
    if (fy > (1.0 - kEdgeThreshold)) {
      has_top = true;
    }
  }

  // If no signal pins were found on any edge, apply halo on all sides
  // as a safe default.
  if (!has_left && !has_right && !has_bottom && !has_top) {
    return {halo_x, halo_y, halo_x, halo_y};
  }

  Halo h;
  h.left = has_left ? halo_x : 0;
  h.right = has_right ? halo_x : 0;
  h.bottom = has_bottom ? halo_y : 0;
  h.top = has_top ? halo_y : 0;
  return h;
}

// Apply orientation transform to a Halo that was computed in R0.
// Macro orientations rotate which physical side corresponds to which
// logical direction.
static Halo orientHalo(const Halo& h, odb::dbOrientType orient)
{
  Halo result;
  switch (orient.getValue()) {
    case odb::dbOrientType::R0:
      result = h;
      break;
    case odb::dbOrientType::R90:
      // 90° CCW: left→bottom, top→left, right→top, bottom→right
      result.left = h.bottom;
      result.bottom = h.right;
      result.right = h.top;
      result.top = h.left;
      break;
    case odb::dbOrientType::R180:
      result.left = h.right;
      result.bottom = h.top;
      result.right = h.left;
      result.top = h.bottom;
      break;
    case odb::dbOrientType::R270:
      result.left = h.top;
      result.bottom = h.left;
      result.right = h.bottom;
      result.top = h.right;
      break;
    case odb::dbOrientType::MY:
      // Mirror about Y-axis: left↔right
      result.left = h.right;
      result.bottom = h.bottom;
      result.right = h.left;
      result.top = h.top;
      break;
    case odb::dbOrientType::MX:
      // Mirror about X-axis: top↔bottom
      result.left = h.left;
      result.bottom = h.top;
      result.right = h.right;
      result.top = h.bottom;
      break;
    default:
      result = h;
      break;
  }
  return result;
}

static Halo computePinAwareHaloForInst(odb::dbInst* inst, int halo_x, int halo_y)
{
  odb::Rect master_box;
  inst->getMaster()->getPlacementBoundary(master_box);
  if (master_box.dx() == 0 || master_box.dy() == 0) {
    Halo h;
    h.left = halo_x;
    h.right = halo_x;
    h.bottom = halo_y;
    h.top = halo_y;
    return h;
  }

  bool has_left = false;
  bool has_right = false;
  bool has_bottom = false;
  bool has_top = false;

  for (odb::dbITerm* iterm : inst->getITerms()) {
    if (iterm->getSigType() == odb::dbSigType::SIGNAL) {
      int x, y;
      odb::Point loc;
      iterm->getAvgXY(&x, &y);
      loc.setX(x);
      loc.setY(y);
      if (loc.x() <= master_box.xMin()) {
        has_left = true;
      }
      if (loc.x() >= master_box.xMax()) {
        has_right = true;
      }
      if (loc.y() <= master_box.yMin()) {
        has_bottom = true;
      }
      if (loc.y() >= master_box.yMax()) {
        has_top = true;
      }
    }
  }

  if (!has_left && !has_right && !has_bottom && !has_top) {
    has_left = has_right = has_bottom = has_top = true;
  }

  Halo h;
  h.left = has_left ? halo_x : 0;
  h.right = has_right ? halo_x : 0;
  h.bottom = has_bottom ? halo_y : 0;
  h.top = has_top ? halo_y : 0;
  return h;
}

void BStarTree::applyDefaultHalo(BStarNode* node)
{
  if (node == nullptr) {
    return;
  }

  if (node->isSoftMacro()) {
    // if (default_halo_configured_) {
    //   Halo h;
    //   h.left = default_halo_x_;
    //   h.right = default_halo_x_;
    //   h.bottom = default_halo_y_;
    //   h.top = default_halo_y_;
    //   node->setHalo(h);
    // }
    return;
  }

  odb::dbInst* inst = node->getInst();
  if (inst == nullptr) {
    return;
  }

  auto it = macro_halos_.find(inst);
  if (it != macro_halos_.end()) {
    Halo oriented = orientHalo(it->second, node->getOrientation());
    node->setHalo(oriented);
    return;
  }

  if (!default_halo_configured_) {
    return;
  }

  if (pin_aware_halo_enabled_) {
    Halo base = computePinAwareHaloForInst(inst, default_halo_x_, default_halo_y_);
    Halo oriented = orientHalo(base, node->getOrientation());
    node->setHalo(oriented);
    return;
  }

  Halo h;
  h.left = default_halo_x_;
  h.right = default_halo_x_;
  h.bottom = default_halo_y_;
  h.top = default_halo_y_;
  node->setHalo(h);
}

void BStarTree::computePinAwareHalos(int halo_x, int halo_y)
{
  if (halo_x <= 0 && halo_y <= 0) {
    return;
  }

  default_halo_x_ = halo_x;
  default_halo_y_ = halo_y;
  default_halo_configured_ = true;
  pin_aware_halo_enabled_ = true;

  for (auto& node : nodes_) {
    if (node->isSoftMacro()) {
      // Soft macros have no halo.
      continue;
    }

    odb::dbInst* inst = node->getInst();

    // Check for per-instance override first
    auto it = macro_halos_.find(inst);
    if (it != macro_halos_.end()) {
      Halo oriented = orientHalo(it->second, node->getOrientation());
      node->setHalo(oriented);
      continue;
    }

    // Compute pin-aware halo in R0, then orient
    Halo base = computePinAwareHaloForInst(inst, halo_x, halo_y);
    Halo oriented = orientHalo(base, node->getOrientation());
    node->setHalo(oriented);
  }
}

void BStarTree::setUniformHalo(int halo_x, int halo_y)
{
  default_halo_x_ = halo_x;
  default_halo_y_ = halo_y;
  default_halo_configured_ = true;
  pin_aware_halo_enabled_ = false;

  for (auto& node : nodes_) {
    if (node->isSoftMacro()) {
      // No per-instance override for soft macros; always use uniform.
      Halo h;
      h.left = halo_x;
      h.right = halo_x;
      h.bottom = halo_y;
      h.top = halo_y;
      node->setHalo(h);
      continue;
    }

    odb::dbInst* inst = node->getInst();

    // Check for per-instance override first
    auto it = macro_halos_.find(inst);
    if (it != macro_halos_.end()) {
      Halo oriented = orientHalo(it->second, node->getOrientation());
      node->setHalo(oriented);
      continue;
    }

    Halo h;
    h.left = halo_x;
    h.right = halo_x;
    h.bottom = halo_y;
    h.top = halo_y;
    node->setHalo(h);
  }
}

void BStarTree::setMacroHalo(odb::dbInst* inst, const Halo& halo)
{
  macro_halos_[inst] = halo;

  // If the node already exists, apply immediately
  for (auto& node : nodes_) {
    if (node->getInst() == inst) {
      Halo oriented = orientHalo(halo, node->getOrientation());
      node->setHalo(oriented);
      break;
    }
  }
}

}  // namespace pne
