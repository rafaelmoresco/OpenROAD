// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

namespace pne {

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

}  // namespace pne
