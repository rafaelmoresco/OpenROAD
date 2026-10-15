// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/ChannelLegalizer.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace pne {

namespace {

int64_t overlapArea(const ChannelBox& a, const ChannelBox& b)
{
  const int64_t dx = static_cast<int64_t>(std::min(a.xhi, b.xhi))
                     - std::max(a.xlo, b.xlo);
  const int64_t dy = static_cast<int64_t>(std::min(a.yhi, b.yhi))
                     - std::max(a.ylo, b.ylo);
  return (dx > 0 && dy > 0) ? dx * dy : 0;
}

int64_t outsideArea(const ChannelBox& box, const ChannelBox& core)
{
  const int64_t area = static_cast<int64_t>(box.xhi - box.xlo)
                       * static_cast<int64_t>(box.yhi - box.ylo);
  return area - overlapArea(box, core);
}

ChannelBox haloBox(const ChannelMacro& macro)
{
  return {macro.body.xlo - macro.halo.left,
          macro.body.ylo - macro.halo.bottom,
          macro.body.xhi + macro.halo.right,
          macro.body.yhi + macro.halo.top};
}

// A third body spanning the whole channel splits it into two channels of
// their own, so the outer pair no longer faces.
bool isSplit(const std::vector<ChannelMacro>& macros,
             const ChannelBox& region,
             bool along_x,
             int skip1,
             int skip2)
{
  for (int k = 0; k < static_cast<int>(macros.size()); ++k) {
    if (k == skip1 || k == skip2) {
      continue;
    }
    const ChannelBox& b = macros[k].body;
    const bool split
        = along_x ? (b.xlo < region.xhi && b.xhi > region.xlo
                     && b.ylo <= region.ylo && b.yhi >= region.yhi)
                  : (b.ylo < region.yhi && b.yhi > region.ylo
                     && b.xlo <= region.xlo && b.xhi >= region.xhi);
    if (split) {
      return true;
    }
  }
  return false;
}

void checkPair(const std::vector<ChannelMacro>& macros,
               int i,
               int j,
               const ChannelRule& rule,
               std::vector<NarrowChannel>& channels)
{
  const ChannelBox& a = macros[i].body;
  const ChannelBox& b = macros[j].body;

  if (rule.min_x > 0) {
    const int ylo = std::max(a.ylo, b.ylo);
    const int yhi = std::min(a.yhi, b.yhi);
    if (yhi > ylo && (a.xhi <= b.xlo || b.xhi <= a.xlo)) {
      const int low = a.xhi <= b.xlo ? i : j;
      const int high = low == i ? j : i;
      const ChannelBox region{
          macros[low].body.xhi, ylo, macros[high].body.xlo, yhi};
      const int gap = region.xhi - region.xlo;
      if (gap > 0 && gap < rule.min_x
          && !isSplit(macros, region, true, i, j)) {
        channels.push_back({low, high, true, gap, rule.min_x, region});
      }
    }
  }

  if (rule.min_y > 0) {
    const int xlo = std::max(a.xlo, b.xlo);
    const int xhi = std::min(a.xhi, b.xhi);
    if (xhi > xlo && (a.yhi <= b.ylo || b.yhi <= a.ylo)) {
      const int low = a.yhi <= b.ylo ? i : j;
      const int high = low == i ? j : i;
      const ChannelBox region{
          xlo, macros[low].body.yhi, xhi, macros[high].body.ylo};
      const int gap = region.yhi - region.ylo;
      if (gap > 0 && gap < rule.min_y
          && !isSplit(macros, region, false, i, j)) {
        channels.push_back({low, high, false, gap, rule.min_y, region});
      }
    }
  }
}

void checkCoreEdges(const std::vector<ChannelMacro>& macros,
                    int i,
                    const ChannelRule& rule,
                    std::vector<NarrowChannel>& channels)
{
  constexpr int kEdge = NarrowChannel::kCoreEdge;
  const ChannelBox& b = macros[i].body;
  const ChannelBox& core = rule.core;

  auto check = [&](bool along_x, int low, int high, const ChannelBox& region) {
    const int gap
        = along_x ? region.xhi - region.xlo : region.yhi - region.ylo;
    const int required = along_x ? rule.min_x : rule.min_y;
    if (required > 0 && gap > 0 && gap < required
        && !isSplit(macros, region, along_x, i, i)) {
      channels.push_back({low, high, along_x, gap, required, region});
    }
  };

  check(true, kEdge, i, {core.xlo, b.ylo, b.xlo, b.yhi});
  check(true, i, kEdge, {b.xhi, b.ylo, core.xhi, b.yhi});
  check(false, kEdge, i, {b.xlo, core.ylo, b.xhi, b.ylo});
  check(false, i, kEdge, {b.xlo, b.yhi, b.xhi, core.yhi});
}

// Deficit of the channel between two bodies; 0 when abutted, wide enough
// or not facing each other.
int pairDeficit(const ChannelBox& a, const ChannelBox& b, const ChannelRule& rule)
{
  if (rule.min_x > 0 && std::min(a.yhi, b.yhi) > std::max(a.ylo, b.ylo)) {
    const int gap = std::max(b.xlo - a.xhi, a.xlo - b.xhi);
    if (gap > 0 && gap < rule.min_x) {
      return rule.min_x - gap;
    }
  }
  if (rule.min_y > 0 && std::min(a.xhi, b.xhi) > std::max(a.xlo, b.xlo)) {
    const int gap = std::max(b.ylo - a.yhi, a.ylo - b.yhi);
    if (gap > 0 && gap < rule.min_y) {
      return rule.min_y - gap;
    }
  }
  return 0;
}

void shiftMacro(ChannelMacro& macro, bool along_x, int delta)
{
  if (along_x) {
    macro.body.xlo += delta;
    macro.body.xhi += delta;
  } else {
    macro.body.ylo += delta;
    macro.body.yhi += delta;
  }
}

// Shifts seed and drags along every macro whose relation to a shifted one
// gets worse, so the group moves rigidly.  Returns the group size, or -1 if
// a fixed macro or obstacle is in the way or a halo would leave the core.
int shiftGroup(std::vector<ChannelMacro>& macros,
               const std::vector<ChannelBox>& obstacles,
               const ChannelRule& rule,
               int seed,
               bool along_x,
               int delta)
{
  const std::vector<ChannelMacro> before = macros;
  std::vector<char> moved(macros.size(), 0);
  std::vector<int> group{seed};
  moved[seed] = 1;
  shiftMacro(macros[seed], along_x, delta);

  for (size_t g = 0; g < group.size(); ++g) {
    const int s = group[g];
    const ChannelBox old_halo = haloBox(before[s]);
    const ChannelBox new_halo = haloBox(macros[s]);
    if (outsideArea(new_halo, rule.core) > outsideArea(old_halo, rule.core)) {
      return -1;
    }
    for (const ChannelBox& obstacle : obstacles) {
      if (overlapArea(new_halo, obstacle) > overlapArea(old_halo, obstacle)) {
        return -1;
      }
    }
    for (int k = 0; k < static_cast<int>(macros.size()); ++k) {
      if (moved[k]) {
        continue;
      }
      const ChannelBox other_halo = haloBox(macros[k]);
      const bool worse
          = overlapArea(new_halo, other_halo) > overlapArea(old_halo, other_halo)
            || pairDeficit(macros[s].body, macros[k].body, rule)
                   > pairDeficit(before[s].body, macros[k].body, rule);
      if (!worse) {
        continue;
      }
      if (!macros[k].movable) {
        return -1;
      }
      moved[k] = 1;
      shiftMacro(macros[k], along_x, delta);
      group.push_back(k);
    }
  }
  return static_cast<int>(group.size());
}

int lo(const ChannelBox& box, bool along_x)
{
  return along_x ? box.xlo : box.ylo;
}

int hi(const ChannelBox& box, bool along_x)
{
  return along_x ? box.xhi : box.yhi;
}

int lowSideHalo(const Halo& halo, bool along_x)
{
  return along_x ? halo.left : halo.bottom;
}

int highSideHalo(const Halo& halo, bool along_x)
{
  return along_x ? halo.right : halo.top;
}

// Distances at which macro m, moving along the axis in direction dir, would
// abut or sit a full channel away from a facing body or the core edge.
std::vector<int> stopDistances(const std::vector<ChannelMacro>& macros,
                               const ChannelRule& rule,
                               int m,
                               bool along_x,
                               int dir)
{
  const ChannelBox& body = macros[m].body;
  const int required = along_x ? rule.min_x : rule.min_y;
  const int lead_halo = dir > 0 ? highSideHalo(macros[m].halo, along_x)
                                : lowSideHalo(macros[m].halo, along_x);
  std::vector<int> stops;
  auto add = [&](int gap, int facing_halo) {
    if (lead_halo == 0 && facing_halo == 0) {
      stops.push_back(gap);
    }
    stops.push_back(gap - required);
  };

  add(dir > 0 ? hi(rule.core, along_x) - hi(body, along_x)
              : lo(body, along_x) - lo(rule.core, along_x),
      0);
  for (int k = 0; k < static_cast<int>(macros.size()); ++k) {
    const ChannelBox& other = macros[k].body;
    if (k == m
        || std::min(hi(other, !along_x), hi(body, !along_x))
               <= std::max(lo(other, !along_x), lo(body, !along_x))) {
      continue;
    }
    if (dir > 0 && lo(other, along_x) >= hi(body, along_x)) {
      add(lo(other, along_x) - hi(body, along_x),
          lowSideHalo(macros[k].halo, along_x));
    } else if (dir < 0 && hi(other, along_x) <= lo(body, along_x)) {
      add(lo(body, along_x) - hi(other, along_x),
          highSideHalo(macros[k].halo, along_x));
    }
  }
  return stops;
}

struct Shift
{
  int macro;
  bool along_x;
  int delta;
};

constexpr int kNumTiers = 3;
constexpr int kMaxStops = 3;

// Ways to fix one channel, cheapest tier first: 0 abuts (if neither side has
// a halo there) or widens it to the minimum; 1 moves a side further, to the
// next abutment or full channel; 2 slides a side across until they stop facing.
std::vector<std::vector<Shift>> fixCandidates(
    const NarrowChannel& channel,
    const std::vector<ChannelMacro>& macros,
    const ChannelRule& rule,
    int tier)
{
  constexpr int kEdge = NarrowChannel::kCoreEdge;
  const bool ax = channel.along_x;
  const int low = channel.low;
  const int high = channel.high;
  const bool low_movable = low != kEdge && macros[low].movable;
  const bool high_movable = high != kEdge && macros[high].movable;
  const int widen = channel.required - channel.gap;

  std::vector<std::vector<Shift>> candidates;
  if (tier == 0) {
    const int low_halo = low == kEdge ? 0 : highSideHalo(macros[low].halo, ax);
    const int high_halo
        = high == kEdge ? 0 : lowSideHalo(macros[high].halo, ax);
    if (low_halo == 0 && high_halo == 0) {
      if (high_movable) {
        candidates.push_back({{high, ax, -channel.gap}});
      }
      if (low_movable) {
        candidates.push_back({{low, ax, channel.gap}});
      }
    }
    if (high_movable) {
      candidates.push_back({{high, ax, widen}});
    }
    if (low_movable) {
      candidates.push_back({{low, ax, -widen}});
    }
    if (low_movable && high_movable) {
      const int grid = std::max(1, rule.grid);
      const int half = (widen / 2) / grid * grid;
      if (half > 0 && half < widen) {
        candidates.push_back({{low, ax, -half}, {high, ax, widen - half}});
      }
    }
  } else if (tier == 1) {
    for (const auto& [m, dir] : {std::pair{high, 1}, std::pair{low, -1}}) {
      if (m == kEdge || !macros[m].movable) {
        continue;
      }
      std::vector<int> stops = stopDistances(macros, rule, m, ax, dir);
      std::sort(stops.begin(), stops.end());
      stops.erase(std::unique(stops.begin(), stops.end()), stops.end());
      int taken = 0;
      for (const int stop : stops) {
        if (stop > widen && taken++ < kMaxStops) {
          candidates.push_back({{m, ax, dir * stop}});
        }
      }
    }
  } else if (low != kEdge && high != kEdge) {
    for (const auto& [m, other] : {std::pair{high, low}, std::pair{low, high}}) {
      if (!macros[m].movable) {
        continue;
      }
      const ChannelBox& body = macros[m].body;
      const ChannelBox& facing = macros[other].body;
      candidates.push_back({{m, !ax, hi(facing, !ax) - lo(body, !ax)}});
      candidates.push_back({{m, !ax, lo(facing, !ax) - hi(body, !ax)}});
    }
  }
  return candidates;
}

// Violations are compared by count first, then by total deficit.
std::pair<size_t, int64_t> severity(const std::vector<NarrowChannel>& channels)
{
  int64_t deficit = 0;
  for (const NarrowChannel& channel : channels) {
    deficit += channel.required - channel.gap;
  }
  return {channels.size(), deficit};
}

}  // namespace

std::vector<NarrowChannel> findNarrowChannels(
    const std::vector<ChannelMacro>& macros,
    const ChannelRule& rule)
{
  std::vector<NarrowChannel> channels;
  if (!rule.enabled()) {
    return channels;
  }

  const int n = static_cast<int>(macros.size());
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&macros](int a, int b) {
    return std::tie(macros[a].body.xlo, a) < std::tie(macros[b].body.xlo, b);
  });

  // Sweep in x: a later body can only face this one through an x overlap or
  // a narrow x gap, so stop once it starts min_x past this body.
  for (int p = 0; p < n; ++p) {
    const int i = order[p];
    const int64_t reach = static_cast<int64_t>(macros[i].body.xhi) + rule.min_x;
    for (int q = p + 1; q < n; ++q) {
      const int j = order[q];
      if (macros[j].body.xlo >= reach) {
        break;
      }
      if (macros[i].movable || macros[j].movable) {
        checkPair(macros, i, j, rule, channels);
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    if (macros[i].movable) {
      checkCoreEdges(macros, i, rule, channels);
    }
  }

  std::sort(channels.begin(),
            channels.end(),
            [](const NarrowChannel& a, const NarrowChannel& b) {
              return std::tie(a.region.xlo, a.region.ylo, a.along_x, a.low, a.high)
                     < std::tie(b.region.xlo, b.region.ylo, b.along_x, b.low, b.high);
            });
  return channels;
}

double narrowChannelPenalty(const std::vector<NarrowChannel>& channels)
{
  double penalty = 0.0;
  for (const NarrowChannel& channel : channels) {
    penalty += static_cast<double>(channel.required - channel.gap)
               / static_cast<double>(channel.required);
  }
  return penalty;
}

std::vector<NarrowChannel> legalizeNarrowChannels(
    std::vector<ChannelMacro>& macros,
    const std::vector<ChannelBox>& obstacles,
    const ChannelRule& rule)
{
  std::vector<NarrowChannel> channels = findNarrowChannels(macros, rule);
  auto current = severity(channels);

  // Every accepted fix strictly lowers the severity, so this terminates; the
  // cap only bounds the runtime on pathological inputs.
  const size_t max_fixes = 4 * channels.size() + 16;
  for (size_t fix = 0; fix < max_fixes && !channels.empty(); ++fix) {
    // Prefer the most local fix anywhere over a bigger move on an earlier
    // channel: escalate a tier only when no channel can be fixed in this one.
    bool accepted = false;
    for (int tier = 0; tier < kNumTiers && !accepted; ++tier) {
      for (const NarrowChannel& channel : channels) {
        bool found = false;
        std::vector<ChannelMacro> best;
        std::vector<NarrowChannel> best_channels;
        auto best_severity = current;
        int64_t best_displacement = 0;

        for (const auto& candidate :
             fixCandidates(channel, macros, rule, tier)) {
          std::vector<ChannelMacro> trial = macros;
          int64_t displacement = 0;
          bool feasible = true;
          for (const Shift& shift : candidate) {
            const int group = shiftGroup(trial,
                                         obstacles,
                                         rule,
                                         shift.macro,
                                         shift.along_x,
                                         shift.delta);
            if (group < 0) {
              feasible = false;
              break;
            }
            displacement
                += static_cast<int64_t>(group) * std::abs(shift.delta);
          }
          if (!feasible) {
            continue;
          }
          std::vector<NarrowChannel> trial_channels
              = findNarrowChannels(trial, rule);
          const auto trial_severity = severity(trial_channels);
          if (!(trial_severity < current)) {
            continue;
          }
          if (!found || trial_severity < best_severity
              || (trial_severity == best_severity
                  && displacement < best_displacement)) {
            found = true;
            best = std::move(trial);
            best_channels = std::move(trial_channels);
            best_severity = trial_severity;
            best_displacement = displacement;
          }
        }

        if (found) {
          macros = std::move(best);
          channels = std::move(best_channels);
          current = best_severity;
          accepted = true;
          break;
        }
      }
    }
    if (!accepted) {
      break;
    }
  }
  return channels;
}

}  // namespace pne
