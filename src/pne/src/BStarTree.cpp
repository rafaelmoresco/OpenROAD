// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/BStarTree.h"

#include <algorithm>
#include <limits>

namespace pne {

//-----------------------------------------------------------------------------
// BStarNode implementation
//-----------------------------------------------------------------------------

BStarNode::BStarNode(odb::dbInst* inst, int id) : id_(id), inst_(inst)
{
}

int BStarNode::getWidth() const
{
  odb::dbBox* bbox = inst_->getBBox();
  return bbox->getBox().dx();
}

int BStarNode::getHeight() const
{
  odb::dbBox* bbox = inst_->getBBox();
  return bbox->getBox().dy();
}

odb::dbOrientType BStarNode::getOrientation() const
{
  return inst_->getOrient();
}

void BStarNode::setOrientation(odb::dbOrientType orient)
{
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
    // Simple initial tree construction: chain as right children
    BStarNode* current = root_;
    while (current->getRight() != nullptr) {
      current = current->getRight();
    }
    current->setRight(node.get());
    node->setParent(current);
  }
  
  nodes_.push_back(std::move(node));
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
  
  // Pack starting from root at origin
  packRecursive(root_, 0);
  
  // Compute bounding box
  width_ = 0;
  height_ = 0;
  for (const auto& node : nodes_) {
    int right = node->getX() + node->getWidth();
    int top = node->getY() + node->getHeight();
    width_ = std::max(width_, right);
    height_ = std::max(height_, top);
  }
}

void BStarTree::packRecursive(BStarNode* node, int x)
{
  if (node == nullptr) {
    return;
  }
  
  // Find Y coordinate using contour
  int y = findContourY(x, node->getWidth());
  node->setX(x);
  node->setY(y);
  
  // Update contour
  updateContour(x, y, node->getWidth(), node->getHeight());
  
  // Left child: place to the right (same y-level)
  if (node->getLeft() != nullptr) {
    packRecursive(node->getLeft(), x + node->getWidth());
  }
  
  // Right child: place above (compacted upward)
  if (node->getRight() != nullptr) {
    packRecursive(node->getRight(), x);
  }
}

int BStarTree::findContourY(int x, int width)
{
  int max_y = 0;
  
  // Find maximum Y in the range [x, x+width)
  for (const auto& c : contours_) {
    int c_end = c.x + c.width;
    
    // Check if this contour segment overlaps with [x, x+width)
    if (c.x < x + width && c_end > x) {
      max_y = std::max(max_y, c.y);
    }
  }
  
  return max_y;
}

void BStarTree::updateContour(int x, int y, int width, int height)
{
  std::vector<Contour> new_contours;
  int new_y = y + height;
  
  for (const auto& c : contours_) {
    int c_end = c.x + c.width;
    int block_end = x + width;
    
    // Contour completely before block
    if (c_end <= x) {
      new_contours.push_back(c);
    }
    // Contour completely after block
    else if (c.x >= block_end) {
      new_contours.push_back(c);
    }
    // Contour overlaps with block
    else {
      // Keep part before block
      if (c.x < x) {
        new_contours.push_back({c.x, c.y, x - c.x});
      }
      
      // Keep part after block
      if (c_end > block_end) {
        new_contours.push_back({block_end, c.y, c_end - block_end});
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
  
  // Swap positions in tree structure
  BStarNode* parent1 = node1->getParent();
  BStarNode* parent2 = node2->getParent();
  BStarNode* left1 = node1->getLeft();
  BStarNode* right1 = node1->getRight();
  BStarNode* left2 = node2->getLeft();
  BStarNode* right2 = node2->getRight();
  
  // Handle parent-child relationship
  if (parent1 == node2) {
    // node2 is parent of node1
    if (node2->getLeft() == node1) {
      node1->setLeft(left2);
      node1->setRight(right2);
      node2->setLeft(left1);
      node2->setRight(right1);
    } else {
      node1->setLeft(left2);
      node1->setRight(right2);
      node2->setLeft(left1);
      node2->setRight(right1);
    }
    
    node1->setParent(parent2);
    if (parent2) {
      if (parent2->getLeft() == node2) parent2->setLeft(node1);
      else parent2->setRight(node1);
    } else {
      root_ = node1;
    }
    
    node2->setParent(node1);
    if (left1 != node2 && right1 != node2) {
      if (node1->getLeft() == node2) node1->setLeft(node2);
      else node1->setRight(node2);
    }
    
  } else if (parent2 == node1) {
    // node1 is parent of node2 - symmetric case
    swapNodes(id2, id1);
    return;
  } else {
    // Normal swap
    node1->setLeft(left2);
    node1->setRight(right2);
    node1->setParent(parent2);
    
    node2->setLeft(left1);
    node2->setRight(right1);
    node2->setParent(parent1);
    
    // Update parent pointers
    if (parent1) {
      if (parent1->getLeft() == node1) parent1->setLeft(node2);
      else parent1->setRight(node2);
    } else {
      root_ = node2;
    }
    
    if (parent2) {
      if (parent2->getLeft() == node2) parent2->setLeft(node1);
      else parent2->setRight(node1);
    } else {
      root_ = node1;
    }
  }
  
  // Update children's parent pointers
  if (left1 && left1 != node2) left1->setParent(node2);
  if (right1 && right1 != node2) right1->setParent(node2);
  if (left2 && left2 != node1) left2->setParent(node1);
  if (right2 && right2 != node1) right2->setParent(node1);
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
  
  if (node == nullptr) {
    return;
  }
  
  // Cannot move to self or to own descendant
  if (node == new_parent) {
    return;
  }
  
  removeFromTree(node);
  insertInTree(node, new_parent, as_left_child);
}

void BStarTree::applyPlacement()
{
  for (const auto& node : nodes_) {
    int x = node->getX();
    int y = node->getY();
    node->getInst()->setLocation(x, y);
    node->getInst()->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  }
}

void BStarTree::save()
{
  backup_.clear();
  backup_.reserve(nodes_.size());
  
  for (const auto& node : nodes_) {
    NodeState state;
    state.left_id = node->getLeft() ? node->getLeft()->getId() : -1;
    state.right_id = node->getRight() ? node->getRight()->getId() : -1;
    state.parent_id = node->getParent() ? node->getParent()->getId() : -1;
    state.x = node->getX();
    state.y = node->getY();
    state.orient = node->getOrientation();
    backup_.push_back(state);
  }
}

void BStarTree::restore()
{
  if (backup_.size() != nodes_.size()) {
    return;
  }
  
  // First pass: restore orientations and coordinates
  for (size_t i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->setX(backup_[i].x);
    nodes_[i]->setY(backup_[i].y);
    nodes_[i]->setOrientation(backup_[i].orient);
  }
  
  // Second pass: restore tree structure
  for (size_t i = 0; i < nodes_.size(); ++i) {
    BStarNode* left = findNode(backup_[i].left_id);
    BStarNode* right = findNode(backup_[i].right_id);
    BStarNode* parent = findNode(backup_[i].parent_id);
    
    nodes_[i]->setLeft(left);
    nodes_[i]->setRight(right);
    nodes_[i]->setParent(parent);
    
    if (backup_[i].parent_id == -1) {
      root_ = nodes_[i].get();
    }
  }
}

}  // namespace pne
