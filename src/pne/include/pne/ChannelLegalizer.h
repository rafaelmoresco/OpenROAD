// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <vector>

#include "pne/Halo.h"

namespace pne {

// Axis-aligned box in DBU.
struct ChannelBox
{
  int xlo = 0;
  int ylo = 0;
  int xhi = 0;
  int yhi = 0;
};

// A hard macro as the channel rule sees it: the body (no halo), the halo
// other macros must stay out of, and whether the macro may be moved.
struct ChannelMacro
{
  ChannelBox body;
  Halo halo;
  bool movable = true;
};

// Facing macro bodies, and a movable body and the core edge, must either
// abut or be at least min_x (gaps along x) / min_y (gaps along y) apart.
struct ChannelRule
{
  ChannelBox core;
  int min_x = 0;  // 0 disables the rule along x
  int min_y = 0;  // 0 disables the rule along y
  int grid = 1;   // legalization shifts stay multiples of this

  bool enabled() const { return min_x > 0 || min_y > 0; }
};

struct NarrowChannel
{
  static constexpr int kCoreEdge = -1;

  int low = kCoreEdge;   // macro on the left/bottom side, or the core edge
  int high = kCoreEdge;  // macro on the right/top side, or the core edge
  bool along_x = true;   // gap measured along x, i.e. a vertical channel
  int gap = 0;
  int required = 0;
  ChannelBox region;  // the sliver between the two sides
};

// Every gap that is neither 0 nor at least the rule minimum, in a fixed
// order.  Pairs of immovable macros are skipped since nothing can fix them.
std::vector<NarrowChannel> findNarrowChannels(
    const std::vector<ChannelMacro>& macros,
    const ChannelRule& rule);

// Sum of (required - gap) / required: one unit per fully closed channel.
double narrowChannelPenalty(const std::vector<NarrowChannel>& channels);

// Abuts or widens narrow channels by rigidly shifting movable macros.  A
// shift is kept only if it adds no halo/obstacle overlap, keeps halos inside
// the core and strictly reduces the violations.  Returns what is left.
std::vector<NarrowChannel> legalizeNarrowChannels(
    std::vector<ChannelMacro>& macros,
    const std::vector<ChannelBox>& obstacles,
    const ChannelRule& rule);

}  // namespace pne
