// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <memory>
#include <vector>

#include "odb/db.h"

namespace pne {

// B*-Tree node representing a macro instance
class BStarNode
{
 public:
  BStarNode(odb::dbInst* inst, int id);

  // Node properties
  int getId() const { return id_; }
  odb::dbInst* getInst() const { return inst_; }
  
  // Tree structure
  BStarNode* getLeft() const { return left_; }
  BStarNode* getRight() const { return right_; }
  BStarNode* getParent() const { return parent_; }
  
  void setLeft(BStarNode* node) { left_ = node; }
  void setRight(BStarNode* node) { right_ = node; }
  void setParent(BStarNode* node) { parent_ = node; }
  
  // Placement coordinates (computed during packing)
  int getX() const { return x_; }
  int getY() const { return y_; }
  void setX(int x) { x_ = x; }
  void setY(int y) { y_ = y; }
  
  // Macro dimensions
  int getWidth() const;
  int getHeight() const;
  
  // Orientation
  odb::dbOrientType getOrientation() const;
  void setOrientation(odb::dbOrientType orient);
  
 private:
  int id_;
  odb::dbInst* inst_;
  
  // B*-Tree structure
  BStarNode* left_ = nullptr;
  BStarNode* right_ = nullptr;
  BStarNode* parent_ = nullptr;
  
  // Placement result
  int x_ = 0;
  int y_ = 0;
};

// B*-Tree structure for macro placement
class BStarTree
{
 public:
  BStarTree();
  ~BStarTree();
  
  // Tree construction
  void addMacro(odb::dbInst* inst);
  void clear();
  
  // Tree operations
  BStarNode* getRoot() const { return root_; }
  const std::vector<std::unique_ptr<BStarNode>>& getNodes() const { return nodes_; }
  int getNumNodes() const { return nodes_.size(); }
  
  // Packing: compute coordinates from B*-Tree structure
  void pack();
  
  // Get bounding box after packing
  int getWidth() const { return width_; }
  int getHeight() const { return height_; }
  int getArea() const { return width_ * height_; }
  
  // Perturbation operators for SA
  void swapNodes(int id1, int id2);
  void rotateNode(int id);
  void moveNode(int id, int new_parent_id, bool as_left_child);
  
  // Apply placement to database
  void applyPlacement();
  
  // Copy/restore operations for SA
  void save();
  void restore();
  
 private:
  BStarNode* root_ = nullptr;
  std::vector<std::unique_ptr<BStarNode>> nodes_;
  
  // Bounding box computed during packing
  int width_ = 0;
  int height_ = 0;
  
  // Contour structure for packing
  struct Contour {
    int x;
    int y;
    int width;
  };
  std::vector<Contour> contours_;
  
  // Backup for restore
  struct NodeState {
    int left_id;
    int right_id;
    int parent_id;
    int x;
    int y;
    odb::dbOrientType orient;
  };
  std::vector<NodeState> backup_;
  
  // Helper methods
  void packRecursive(BStarNode* node, int x);
  int findContourY(int x, int width);
  void updateContour(int x, int y, int width, int height);
  BStarNode* findNode(int id);
  void removeFromTree(BStarNode* node);
  void insertInTree(BStarNode* node, BStarNode* parent, bool as_left);
};

}  // namespace pne
