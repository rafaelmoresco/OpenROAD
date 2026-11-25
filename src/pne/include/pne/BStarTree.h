// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "odb/db.h"
#include "pne/Halo.h"

namespace pne {

struct SoftMacro;

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
  int getBodyX() const { return x_ + halo_.left; }
  int getBodyY() const { return y_ + halo_.bottom; }
  
  // Macro dimensions (without halo)
  int getMacroWidth() const;
  int getMacroHeight() const;

  // Effective dimensions (macro + halo)
  int getWidth() const;
  int getHeight() const;
  
  // Orientation.  Held locally for BOTH hard and soft macros: SA flip moves
  // must not touch the DB mid-anneal (parallel multi-start chains share the
  // dbInst objects read-only).  The DB orientation is written only by
  // applyPlacement().
  odb::dbOrientType getOrientation() const;
  void setOrientation(odb::dbOrientType orient);

  // Halo
  void setHalo(const Halo& halo) { halo_ = halo; }
  const Halo& getHalo() const { return halo_; }

  // Index of the subtree (BStarTree root) the node belongs to.
  int getSubtree() const { return subtree_; }
  void setSubtree(int subtree) { subtree_ = subtree; }
  
 private:
  int id_;
  odb::dbInst* inst_;
  SoftMacro* soft_macro_ = nullptr;
  int subtree_ = 0;

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

  static constexpr int kNumCorners = 4;

  BStarTree();
  ~BStarTree();
  
  // Tree construction.
  // Append soft macros using a structured initial placement strategy.
  // This keeps the soft macro subtree within the starting placement region
  // instead of appending them with heap-based insertion.  In a forest they
  // form the soft subtree instead.
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
  // A forest instead (README: set_pine_mp_corner_subtrees): subtree c holds
  // corner_macros[c], built as above within the region and packed toward
  // corner c; an empty soft subtree, centered on the core, follows them.
  void buildForest(
      const std::array<std::vector<odb::dbInst*>, kNumCorners>& corner_macros,
      int region_width,
      int region_height,
      int halo_x = 0,
      int halo_y = 0);
  void clear();
  
  // Tree operations
  BStarNode* getRoot() const { return subtrees_[0].root; }
  const std::vector<std::unique_ptr<BStarNode>>& getNodes() const { return nodes_; }
  int getNumNodes() const { return nodes_.size(); }
  
  // Packing: compute coordinates from B*-Tree structure
  void pack();

  // Subtrees: one for the classic tree, kNumCorners + 1 for a forest.
  bool isForest() const { return subtrees_.size() > 1; }
  int getNumSubtrees() const { return static_cast<int>(subtrees_.size()); }
  // The soft subtree of a forest, or -1.
  int getSoftSubtree() const { return soft_subtree_; }
  BStarNode* getSubtreeRoot(int subtree) const
  {
    return subtrees_[subtree].root;
  }
  // The corner the subtree is packed toward once the layout is anchored.
  Anchor getSubtreeCorner(int subtree) const;
  bool isSubtreeCentered(int subtree) const
  {
    return subtrees_[subtree].centered;
  }
  // Footprint of the subtree's packing (anchor-invariant).
  int getSubtreeWidth(int subtree) const { return subtrees_[subtree].width; }
  int getSubtreeHeight(int subtree) const
  {
    return subtrees_[subtree].height;
  }

  // Corner anchoring configuration.  The core dimensions are needed to
  // reflect the packed layout toward the right / top edges; set them once
  // the placement core is known.  setAnchor selects which corner pack()
  // compacts toward (default BOTTOM_LEFT reproduces the classic packing);
  // in a forest it reflects the whole layout.
  void setCoreDimensions(int width, int height)
  {
    core_width_ = width;
    core_height_ = height;
  }
  void setAnchor(Anchor anchor) { anchor_ = anchor; }
  Anchor getAnchor() const { return anchor_; }
  int getCoreWidth() const { return core_width_; }
  int getCoreHeight() const { return core_height_; }

  // Minimum channel between facing hard macro bodies (DBU, 0 = off).  pack()
  // keeps vertical gaps and gaps to the anchored core edges 0 or this wide.
  void setMinChannel(int min_x, int min_y)
  {
    min_channel_x_ = min_x;
    min_channel_y_ = min_y;
  }

  // Std-cell row grid in the tree frame (y = 0 is the core bottom edge).
  // pack() raises every hard macro until its body, after anchoring, starts
  // at origin + k * pitch.  A pitch of 0 turns it off.
  void setRowGrid(int pitch, int origin)
  {
    row_pitch_ = pitch;
    row_origin_ = origin;
  }

  // Get bounding box after packing.  This is always the true cluster
  // footprint (invariant under the anchor), not the anchored extent.  In a
  // forest it is the extent of the whole layout.
  int getWidth() const { return width_; }
  int getHeight() const { return height_; }
  int64_t getArea() const
  {
    return static_cast<int64_t>(width_) * static_cast<int64_t>(height_);
  }
  
  // Perturbation operators for SA
  void swapNodes(int id1, int id2);
  void rotateNode(int id);
  // The node joins the new parent's subtree.
  void moveNode(int id, int new_parent_id, bool as_left_child);
  // Make the node the subtree's root; the old root becomes its left child.
  void moveToSubtree(int id, int subtree);
  
  // Apply placement to database (writes macro position without halo)
  void applyPlacement();

  // Compute and assign per-node halos based on pin locations.
  // Only sides with signal pins receive the halo margin.
  void computePinAwareHalos(int halo_x, int halo_y);

  // Set a uniform halo on every node (all four sides).
  void setUniformHalo(int halo_x, int halo_y);

  // Set a per-instance halo (all four sides).
  void setMacroHalo(odb::dbInst* inst, const Halo& halo);
  
  // Copy/restore operations for SA.
  void saveSnapshot(SnapshotSlot slot);
  void restoreSnapshot(SnapshotSlot slot);

  // Deep copy for parallel multi-start SA chains: same insts/soft macros
  // (shared, read-only during SA), own structure, coordinates, orientations,
  // halos, and configuration.  Snapshots are not copied.
  std::unique_ptr<BStarTree> clone() const;
  // Adopt the placement state (structure, orientations, halos, anchor) of a
  // clone with identical node ids; the caller re-packs afterwards.
  void copyStateFrom(const BStarTree& other);
  
 private:
  struct Subtree {
    BStarNode* root = nullptr;
    Anchor corner = Anchor::BOTTOM_LEFT;  // before the whole-layout anchor
    bool centered = false;  // centered on the core instead of a corner
    int width = 0;          // footprint
    int height = 0;
  };
  std::vector<Subtree> subtrees_ = std::vector<Subtree>(1);
  int soft_subtree_ = -1;
  std::vector<std::unique_ptr<BStarNode>> nodes_;
  
  // Bounding box computed during packing
  int width_ = 0;
  int height_ = 0;
  
  // Contour structure for packing: the skyline of the halo boxes, and below it
  // that of the bodies.  A body must clear the first and a halo the second, so
  // the halos of neighbours can overlap.
  struct Contour {
    int x;
    int y;  // top of the halo boxes
    int width;
    int body_y = 0;        // top of the bodies (or the floor)
    bool channel = false;  // body_y is the floor or a hard macro
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
    int subtree;
  };
  std::array<std::vector<NodeState>, 3> backups_;
  struct SubtreeState {
    int root_id;
    int width;
    int height;
  };
  std::array<std::vector<SubtreeState>, 3> backup_subtrees_;
  std::array<std::pair<int, int>, 3> backup_sizes_{};
  // The anchor is a tree-level property, so it is snapshotted separately
  // from the per-node state: restoring must reproduce both the arrangement
  // and the corner it was compacted toward.
  std::array<Anchor, 3> backup_anchors_
      = {Anchor::BOTTOM_LEFT, Anchor::BOTTOM_LEFT, Anchor::BOTTOM_LEFT};

  // Corner anchoring state
  Anchor anchor_ = Anchor::BOTTOM_LEFT;
  int core_width_ = 0;
  int core_height_ = 0;

  int min_channel_x_ = 0;
  int min_channel_y_ = 0;
  bool enforce_channels_ = false;  // set per packContour() pass
  // Reflection of the subtree being packed, set per subtree by pack().
  bool flip_x_ = false;
  bool flip_y_ = false;
  int row_pitch_ = 0;
  int row_origin_ = 0;

  // Per-instance halo overrides (set via TCL before placement)
  std::unordered_map<odb::dbInst*, Halo> macro_halos_;
  int default_halo_x_ = 0;
  int default_halo_y_ = 0;
  bool default_halo_configured_ = false;
  bool pin_aware_halo_enabled_ = false;
  
  // Helper methods
  // Tall macros as a left-child chain, the rest in right-child columns, as
  // the nodes of the subtree.
  void buildColumns(const std::vector<odb::dbInst*>& macros,
                    int die_width,
                    int die_height,
                    int halo_x,
                    int halo_y,
                    int subtree);
  void packSubtree(int subtree);
  // Move the bottom-left packing to the subtree's place: reflected toward its
  // corner, or centered on the core.  Leaves the footprint untouched.
  void applyAnchor(int subtree);
  // Bottom-left contour packing and footprint of the subtree, with or
  // without the channel spacing (row alignment applies either way).
  void packContour(int subtree, bool enforce_channels);
  bool flipX() const;
  bool flipY() const;
  // The halo as it faces neighbours in the packed frame: reflection moves the
  // halo box but the body keeps its offset inside it, so the sides swap.
  Halo packedHalo(const BStarNode* node) const;
  // Body x span in the packed frame.
  void bodySpan(const BStarNode* node, int x, int& lo, int& hi) const;
  int edgeChannelOffset(const BStarNode* node) const;
  // Packed x of a node placed right of left: the bodies keep the larger of the
  // facing halos apart, or a full channel.
  int besideX(const BStarNode* left, const BStarNode* node) const;
  // Lowest packed y >= y that puts the node's body on the row grid once
  // anchored.
  int64_t alignToRow(const BStarNode* node, int64_t y) const;
  // left is the node this one is placed right of, or null at the core edge.
  void packRecursive(BStarNode* node,
                     const BStarNode* left,
                     std::vector<bool>& visited);
  int findContourY(int x, int width, const BStarNode* node);
  void updateContour(int x, int y, int width, int height, const BStarNode* node);
  BStarNode* findNode(int id);
  // Take the node out of the tree, its left child (or right) taking its slot.
  void detachNode(BStarNode* node);
  void removeFromTree(BStarNode* node);
  void insertInTree(BStarNode* node, BStarNode* parent, bool as_left);
  void applyDefaultHalo(BStarNode* node);
};

}  // namespace pne
