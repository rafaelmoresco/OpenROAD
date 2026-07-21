// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "odb/db.h"

namespace pne {

struct SoftMacro;

// Per-side halo margins around a macro (in DBU).
// A halo reserves buffer space for routing/timing repair.
struct Halo
{
  int left = 0;
  int bottom = 0;
  int right = 0;
  int top = 0;

  bool hasNonZero() const
  {
    return left > 0 || bottom > 0 || right > 0 || top > 0;
  }
};

// B*-Tree node representing either a hard macro (dbInst) or a soft macro
// (SoftMacro). Hard macro nodes carry a non-null inst_; soft macro nodes
// carry a non-null soft_macro_.
class BStarNode
{
 public:
  // Hard-macro constructor.
  BStarNode(odb::dbInst* inst, int id);
  // Soft-macro constructor.  inst_ is kept null.
  BStarNode(SoftMacro* sm, int id);

  // Node properties
  int getId() const { return id_; }
  odb::dbInst* getInst() const { return inst_; }
  SoftMacro* getSoftMacro() const { return soft_macro_; }

  bool isHardMacro() const { return inst_ != nullptr; }
  bool isSoftMacro() const { return soft_macro_ != nullptr; }
  
  // Tree structure
  BStarNode* getLeft() const { return left_; }
  BStarNode* getRight() const { return right_; }
  BStarNode* getParent() const { return parent_; }
  
  void setLeft(BStarNode* node) { left_ = node; }
  void setRight(BStarNode* node) { right_ = node; }
  void setParent(BStarNode* node) { parent_ = node; }
  
  // Placement coordinates (computed during packing).
  // These include the halo offset: the actual macro origin is at
  // (getX() + halo.left, getY() + halo.bottom).
  int getX() const { return x_; }
  int getY() const { return y_; }
  void setX(int x) { x_ = x; }
  void setY(int y) { y_ = y; }
  
  // Macro dimensions (without halo)
  int getMacroWidth() const;
  int getMacroHeight() const;

  // Effective dimensions (macro + halo)
  int getWidth() const;
  int getHeight() const;
  
  // Orientation
  odb::dbOrientType getOrientation() const;
  void setOrientation(odb::dbOrientType orient);

  // Halo
  void setHalo(const Halo& halo) { halo_ = halo; }
  const Halo& getHalo() const { return halo_; }
  
 private:
  int id_;
  odb::dbInst* inst_;
  SoftMacro* soft_macro_ = nullptr;

  // Cached orientation for soft macro nodes (hard macros read/write the DB).
  odb::dbOrientType local_orient_;
  
  // B*-Tree structure
  BStarNode* left_ = nullptr;
  BStarNode* right_ = nullptr;
  BStarNode* parent_ = nullptr;
  
  // Placement result
  int x_ = 0;
  int y_ = 0;

  // Per-node halo (oriented: accounts for current macro orientation)
  Halo halo_;
};

// B*-Tree structure for macro placement
class BStarTree
{
 public:
  enum class SnapshotSlot {
    CURRENT = 0,
    BEST = 1,
    GLOBAL = 2
  };

  // Corner the packed layout is compacted toward within the core outline.
  // A B*-tree naturally compacts to the bottom-left; reflecting that packing
  // horizontally / vertically yields the equivalent compaction from the other
  // three corners (Chen & Chang, "Modern floorplanning based on B*-tree and
  // fast simulated annealing").  The four candidates differ in how the macro
  // cluster aligns with the fixed IO pins, so they carry different
  // wirelength while sharing the same overlap-free arrangement and footprint.
  enum class Anchor {
    BOTTOM_LEFT = 0,
    BOTTOM_RIGHT = 1,
    TOP_LEFT = 2,
    TOP_RIGHT = 3
  };

  BStarTree();
  ~BStarTree();
  
  // Tree construction
  void addMacro(odb::dbInst* inst);
  // Add a soft macro node using heap-based balanced insertion.
  // The SoftMacro object must outlive the tree (owned by SoftMacroMgr).
  void addSoftMacro(SoftMacro* sm);
  // Append soft macros using a structured initial placement strategy.
  // This keeps the soft macro subtree within the starting placement region
  // instead of appending them with heap-based insertion.
  void addSoftMacros(const std::vector<SoftMacro*>& soft_macros,
                     int die_width,
                     int die_height,
                     int halo_x = 0,
                     int halo_y = 0);
  // Smart initial tree: tall macros as horizontal chain, short macros in columns.
  // halo_y is the per-side vertical halo (top and bottom) that will be applied
  // after tree construction; accounting for it here prevents the initial packing
  // from already exceeding the placement boundary before SA even starts.
  // halo_x is included so horizontal placement capacity estimates honor the
  // available width.
  void buildFromMacros(const std::vector<odb::dbInst*>& macros,
                       int die_width,
                       int die_height,
                       int halo_x = 0,
                       int halo_y = 0);
  void clear();
  
  // Tree operations
  BStarNode* getRoot() const { return root_; }
  const std::vector<std::unique_ptr<BStarNode>>& getNodes() const { return nodes_; }
  int getNumNodes() const { return nodes_.size(); }
  
  // Packing: compute coordinates from B*-Tree structure
  void pack();

  // Corner anchoring configuration.  The core dimensions are needed to
  // reflect the packed layout toward the right / top edges; set them once
  // the placement core is known.  setAnchor selects which corner pack()
  // compacts toward (default BOTTOM_LEFT reproduces the classic packing).
  void setCoreDimensions(int width, int height)
  {
    core_width_ = width;
    core_height_ = height;
  }
  void setAnchor(Anchor anchor) { anchor_ = anchor; }
  Anchor getAnchor() const { return anchor_; }

  // Get bounding box after packing.  This is always the true cluster
  // footprint (invariant under the anchor), not the anchored extent.
  int getWidth() const { return width_; }
  int getHeight() const { return height_; }
  int64_t getArea() const
  {
    return static_cast<int64_t>(width_) * static_cast<int64_t>(height_);
  }
  
  // Perturbation operators for SA
  void swapNodes(int id1, int id2);
  void rotateNode(int id);
  void moveNode(int id, int new_parent_id, bool as_left_child);
  
  // Apply placement to database (writes macro position without halo)
  void applyPlacement();

  // Compute and assign per-node halos based on pin locations.
  // Only sides with signal pins receive the halo margin.
  void computePinAwareHalos(int halo_x, int halo_y);

  // Set a uniform halo on every node (all four sides).
  void setUniformHalo(int halo_x, int halo_y);

  // Set a per-instance halo (all four sides).
  void setMacroHalo(odb::dbInst* inst, const Halo& halo);
  
  // Copy/restore operations for SA (defaults to CURRENT slot).
  void save();
  void restore();
  void saveSnapshot(SnapshotSlot slot);
  void restoreSnapshot(SnapshotSlot slot);
  
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
    Halo halo;
  };
  std::array<std::vector<NodeState>, 3> backups_;
  // The anchor is a tree-level property, so it is snapshotted separately
  // from the per-node state: restoring must reproduce both the arrangement
  // and the corner it was compacted toward.
  std::array<Anchor, 3> backup_anchors_
      = {Anchor::BOTTOM_LEFT, Anchor::BOTTOM_LEFT, Anchor::BOTTOM_LEFT};

  // Corner anchoring state
  Anchor anchor_ = Anchor::BOTTOM_LEFT;
  int core_width_ = 0;
  int core_height_ = 0;

  // Per-instance halo overrides (set via TCL before placement)
  std::unordered_map<odb::dbInst*, Halo> macro_halos_;
  int default_halo_x_ = 0;
  int default_halo_y_ = 0;
  bool default_halo_configured_ = false;
  bool pin_aware_halo_enabled_ = false;
  
  // Helper methods
  // Reflect the bottom-left packing toward the configured corner.  Leaves
  // width_/height_ (the true cluster footprint) untouched.
  void applyAnchor();
  void packRecursive(BStarNode* node, int x, std::vector<bool>& visited);
  int findContourY(int x, int width);
  void updateContour(int x, int y, int width, int height);
  BStarNode* findNode(int id);
  void removeFromTree(BStarNode* node);
  void insertInTree(BStarNode* node, BStarNode* parent, bool as_left);
  void applyDefaultHalo(BStarNode* node);
};

}  // namespace pne
