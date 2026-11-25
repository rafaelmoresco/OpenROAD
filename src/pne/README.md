# PineMP - Pin-Aware Iterative Macro Placer

PineMP (Pin-aware Iterative Macro Placer) is a research macro placement engine that integrates pin assignment inside the macro placer itself and iteratively co-optimizes macro placement and pin assignment.

## Research Goals

Unlike traditional flows that treat pin assignment and macro placement as separate steps, PineMP explores the interaction between these two problems through iterative co-optimization. The key research idea is to adaptively rebalance the cost function between:

- **Internal macro connectivity** - nets connecting macros to each other
- **IO connectivity** - nets connecting macros to IO pads

## Algorithm Overview

1. **Initial Setup**: Random/default pin assignment for all macros
2. **Iterative Co-optimization Loop**:
   - Run simulated annealing (SA) macro placement using B*-Tree representation
   - Cost function uses weighted HPWL: `α * Internal_WL + β * IO_WL + penalties` (or, with the [spread objective](#set_pine_mp_objective), the soft macro part of each net minus a reward for the hard macro part)
   - Initially α (internal weight) is high, β (IO weight) is low
   - After SA convergence: reassign pins based on current placement, with `ppl` seeing each standard cell at its soft macro's center ([seeding](#standard-cell-seeding))
   - Gradually increase β and decrease α
   - Repeat until weights reach target balance (e.g., 50/50)

3. **Final Output**: Optimized macro placement with pin assignments

## Commands

### pine_mp
Performs macro placement using the PineMP algorithm.

```tcl
pine_mp [-num_threads num_threads]
```

**Example:**
```tcl
pine_mp -num_threads 4
```

With `-num_threads N` (N > 1), every outer iteration runs N **parallel multi-start SA chains** — independently seeded annealers on cloned B*-trees — and adopts the lowest-cost result (`PNE-0132` reports the winning chain and the best/worst spread). Single-run SA results vary by up to ~14% between seeds on the benchmark set, so multi-start converts that variance into systematic gains at roughly constant wall-clock time (the chains run concurrently; only chain 0 prints SA progress). With `-num_threads 1` (default) the legacy single-chain behavior is preserved exactly. Chain seeds are deterministic per (iteration, chain), so runs are reproducible.

### set_pine_mp_iterations
Set the number of co-optimization iterations (default: 5).

```tcl
set_pine_mp_iterations num_iterations
```

### set_pine_mp_initial_weights
Set the initial weights for internal and IO nets.

```tcl
set_pine_mp_initial_weights -internal_weight <weight> -io_weight <weight>
```

### set_pine_mp_final_weights
Set the final weights for internal and IO nets.

```tcl
set_pine_mp_final_weights -internal_weight <weight> -io_weight <weight>
```

### set_pine_mp_objective
Choose the wirelength objective (default: `min_wl`).

```tcl
set_pine_mp_objective -style min_wl
set_pine_mp_objective -style spread [-spread_weight weight]
```

`min_wl` minimizes the HPWL of every net. That pulls the hard macros toward the logic they connect to, into the middle of the core, where they block routing and cause congestion. `spread` keeps the standard-cell clusters (soft macros) together and close to the IOs, but pushes the hard macros away from them and from the IOs. It follows TritonMP's `corner_max_wl` style, which maximized macro wirelength.

Negating the HPWL of every net would push the soft macros apart as well, so `spread` splits each net instead. A net's pins are the hard macro pins H, the soft macro centers S and the fixed points F (IO pins, pads, placed cells):

- The **soft part**, `HPWL(S ∪ F)`, is minimized with the internal and IO weights, as `min_wl` minimizes the whole net.
- The **hard extension**, `HPWL(S ∪ F ∪ H) − HPWL(S ∪ F)`, is how far the hard pins stretch the net beyond its soft part. It is rewarded with weight γ (`-spread_weight`, default 1.0).

```
cost = w_int * Σ soft_int / B_int + w_io * Σ soft_io / B_io
     − γ * Σ extension / (N_hard * (W_core + H_core)) + penalties
```

`B_int` and `B_io` are the soft parts of the initial placement. `N_hard` counts the nets with a hard pin and at least one other point, so the reward is measured in core half-perimeters per net. It is at most γ while every net's soft part lies inside the core. Hard pins are clamped to the core first, so a macro gains nothing by overflowing it. The cost can go negative, which the annealer handles, since it only compares costs.

How the pins move:

- **Hard pins** change only the extension, so they are always pushed away from their nets' soft macros, IOs and other hard pins. A hard pin inside its net's soft part feels no push. This matters little for the two- and three-pin nets most macros have.
- **Soft pins** are never pushed out of their net's soft part. Moving one that sets an edge outward always costs, and costs more when a hard pin lies beyond it.
- **Side effect:** on a net with no fixed point, the soft part can move as a whole without growing, and the reward then pushes it away from the hard pins too. A soft macro connected only to a memory drifts away from it.

`PNE-0043` and the final results still report the whole HPWL, which `spread` is expected to raise. `PNE-0048` logs γ and `N_hard`. `PNE-0049` (each iteration) and `PNE-0059` (final results) log the soft parts and the hard extension. [Adaptive weighting](#set_pine_mp_adaptive_weighting) balances the soft parts.

**Limitation:** the B*-tree packs every macro into one cluster compacted toward one corner. So `spread` alone cannot place the hard macros in the four corners around a central block of soft macros. [Corner subtrees](#set_pine_mp_corner_subtrees) lift that limit.

### set_pine_mp_corner_subtrees
Pack the hard macros in four corner subtrees around a centered subtree of soft macros (default: disabled).

```tcl
set_pine_mp_corner_subtrees -enable [-soft_overlap_weight weight]
set_pine_mp_corner_subtrees -disable
```

It is meant for [`-style spread`](#set_pine_mp_objective). A single B*-tree compacts every node toward the bottom-left and then reflects the whole cluster into one corner, so the hard macros cannot reach all four corners and a soft macro reaches the middle only if other nodes fill the space below and left of it. With corner subtrees, the tree has five roots instead of one:

- **Four corner subtrees** hold the hard macros. Each is packed like the single tree, with its own contour, and reflected into its own corner of the core. Shared halos, the channel spacing and row alignment hold within each subtree exactly as before.
- **One centered subtree** holds the soft macros. It is packed the same way and then shifted so that its bounding box is centered on the core.

Packing keeps the macros of one subtree apart but not those of different subtrees; the cost does that:

- **Outline:** each subtree's footprint is checked against the core on its own, and the outline penalty is the sum.
- **Overlap:** the overlap penalty already covers every pair of macros. A soft macro over a hard macro counts only `-soft_overlap_weight` (default 0.1) times its area. Soft macros are only seeds for `ppl` and `gpl`, so the centered block may overlap the corners a little rather than distort them. `0` ignores such overlaps, `1` counts them like any other.
- **Channels:** the [channel rule](#minimum-channel-width) and its legalizer already cover every pair of hard macros.

Hard macros still overlapping in the final placement get a `PNE-0135` warning.

**Initial state.** Largest first, each hard macro goes to the corner farthest from the centroid of the IO pins and pads on its nets, unless it would take that corner past a quarter of the total macro area (an empty corner always takes it). Only pins placed before `pine_mp` count; macros without any go to the least loaded corner. `PNE-0133` logs how many macros each corner got. Each corner builds its columns as the single tree does, in a quarter of the core, and the soft subtree starts as a block with the core's aspect ratio.

**Annealing.**

- Swaps and moves pair hard macros with hard macros, in any corner, and soft macros with soft macros, so each kind stays in its subtrees.
- One move in ten sends a random hard macro to another corner: under a random macro there, or as the root of an empty corner.
- [Slack moves](#set_pine_mp_slack_moves) are off, because they assume one cluster.
- [Corner anchoring](#set_pine_mp_anchoring) still tries four reflections, now of the whole layout, so they swap the corner subtrees with each other.

**Reporting.** `PNE-0134` logs each subtree's corner, number of macros and footprint after every iteration and with the final results, for example `Subtrees: bottom-left 3 (120.0 x 80.5 um), ..., center 12 (310.0 x 280.0 um)`. `PNE-0044` and the final floorplan width, height and area then measure the bounding box of all macros, which spans the core once opposite corners are occupied.

**Open choices:**

- The soft block is centered on the core, not on a point chosen against the IO pins.
- The hard macro pins of a net still push each other apart under `spread`, so macros that talk to each other may end up in different corners.

### set_pine_mp_ppl_macro_hiding
Hide the hard macros from `ppl` during `pine_mp` (default: disabled).

```tcl
set_pine_mp_ppl_macro_hiding -enable
set_pine_mp_ppl_macro_hiding -disable
```

`ppl` places each IO pin close to the instances on its net, hard macros included, while [`-style spread`](#set_pine_mp_objective) pushes the hard macros away from their IO pins, so across iterations the two chase each other. With hiding enabled, PineMP marks the hard macros `UNPLACED` for every `ppl` call and restores their status afterwards. `PNE-0065` logs how many it hid.

`ppl` reads an unplaced instance at the die center, so a hidden macro is not ignored but moved there. A pin whose net reaches only hard macros goes to the middle of the die edge nearest the center. Macros touching the die edge stay visible, because `ppl` keeps pins off the edge along them; this only happens when the core reaches the die edge. The status changes fire the DB callbacks, which only the GUI and the RUDY heat map listen to.

### set_pine_mp_sa_params
Configure simulated annealing parameters.

```tcl
set_pine_mp_sa_params [-initial_temp temp] [-cooling_rate rate] [-max_iterations iterations]
```

The `-initial_temp` and `-cooling_rate` apply to the geometric cooling schedule (used when Fast-SA is disabled); `-max_iterations` bounds the total moves per SA run under either schedule.

### set_pine_mp_slack_moves
Enable or disable slack-based move selection (default: enabled).

```tcl
set_pine_mp_slack_moves -enable [-probability p]
set_pine_mp_slack_moves -disable
```

Following Adya & Markov (ParquetFP, "Fixed-outline floorplanning: enabling hierarchical design"), each block is assigned a *slack* in x and y — how far it can shift toward the right / top before hitting a neighbour or the used extent. Zero-slack blocks are *critical*: they set the floorplan's width / height. A fraction `p` of SA moves (default 0.5) are biased to pick a critical block in the violating dimension and relocate it toward a block with spare room (high slack); the rest stay uniform-random for ergodicity.

The bias is a **feasibility-recovery mechanism**, gated on the *hard*-macro outline with a latch: it engages only while the run has not yet found a placement whose hard macros fit inside the core, and it is **permanently released for the rest of the run** as soon as any best-so-far solution is hard-feasible (`PNE-0038` reports the release point). A rescue budget bounds the attempt: if feasibility has not been reached by 20% of `max_iterations`, the bias gives up and releases anyway (`PNE-0039`), with a fresh no-improvement budget for the uniform phase — this guards against the observed deadlock where the bias itself prevents reaching the fit (all-biased runs on swerv_wrapper froze at an infeasible arrangement that uniform moves resolve). The latch is evaluated on the best-so-far state, not the instantaneous one — the current state wanders through overflowing arrangements throughout the anneal, so an instantaneous gate stays stuck on. After release, moves revert to uniform and the outline penalty alone maintains feasibility, so wirelength optimization is not disturbed. Soft macros are virtual (no DB blockage, no PDN constraint) and are excluded from the gate, though they still act as obstacles in the slack computation.

This staging matters, per benchmarking: an always-on bias fixes designs whose macros overflow (bp: an out-of-bounds/PDN failure becomes a valid fit with lower wirelength) but freezes the search and hurts wirelength on designs that fit easily (bp_multi, swerv); an instantaneous gate reproduced the always-on behavior wherever transient overflow persisted. The best observed behavior is the hybrid: slack packs while infeasible, then uniform moves optimize wirelength (on bp this beat both extremes). Slack is recomputed once per temperature step. Set `-probability 0.0` (or `-disable`) to recover pure uniform-random moves.

Slack moves are always off with [corner subtrees](#set_pine_mp_corner_subtrees), since slack assumes a single packed cluster.

### set_pine_mp_fast_sa
Enable or disable the Fast-SA three-stage annealing schedule (default: **disabled**, experimental).

```tcl
set_pine_mp_fast_sa -enable [-accept_prob P] [-c c] [-k k]
set_pine_mp_fast_sa -disable
```

Fast-SA (Chen & Chang, "Modern floorplanning based on B*-tree and fast simulated annealing") replaces geometric cooling with an adaptive three-stage schedule. With `n` the temperature-step index (each step is a batch of moves) and `delta_cost` the average cost change measured over the previous step:

- **n = 1** — `T1 = delta_avg / -ln(P)`: a high temperature (P near 1) at which almost every uphill move is accepted (random search).
- **2 ≤ n ≤ k** — `Tn = T1 * delta_cost / (n*c)`: the large `c` drives the temperature toward zero, a pseudo-greedy dive to a local minimum.
- **n > k** — `Tn = T1 * delta_cost / n`: dropping `c` makes the temperature jump back up (re-heat) to escape the local minimum, then decay as 1/n while hill-climbing.

**Known limitation (why it is disabled by default):** the literal schedule is mis-scaled for PineMP's penalty-heavy normalized cost. `T1 = delta_avg / -ln(0.99)` is ~100× `delta_avg`, and the stage-3 `1/n` decay is far too slow to reach the low temperatures the geometric schedule uses to pack tightly. Combined with stiff overlap/outline penalties that keep `delta_cost` large, the search never cools: it random-walks at high temperature and the macro cluster overflows the core outline (observed on bp_multi: ~21% height overflow, PDN failure). The geometric auto-calibrated schedule is the working default. Fast-SA is retained for research; making it competitive requires a faster stage-3 cooling (e.g. geometric within stage 3) and/or softer penalties so `delta_cost` shrinks near convergence. Parameters: `-accept_prob` (P, default 0.99), `-c` (default 100), `-k` (default 7).

### set_pine_mp_pin_strategy
Set the pin assignment strategy (uniform, connectivity, random, hungarian).

```tcl
set_pine_mp_pin_strategy connectivity
```

### set_pine_mp_halo
Configure the halo padding around macros. By default, halo is pin-aware and only applies to sides with pins.

```tcl
set_pine_mp_halo -halo_x 10 -halo_y 10
```

Use `-uniform` to force the same halo on all sides regardless of pin location:

```tcl
set_pine_mp_halo -halo_x 10 -halo_y 10 -uniform
```

The halos of neighbouring macros may overlap. No macro body may enter another macro's halo, so facing macros sit the larger of their two facing halos apart rather than the sum. See [Shared halos](#shared-halos), including how this interacts with the minimum channel.

`-min_channel {x y}` (microns; one value applies to both axes) sets the narrowest gap allowed between facing hard macros and between a macro and the core edge, unless they abut. It defaults to twice the halo on each axis (`{20 20}` above), and `0` turns the rule off for that axis. While [row alignment](#set_pine_mp_row_alignment) is on, the y value is rounded up to whole rows (20 µm becomes 21 µm on nangate45's 1.4 µm rows):

```tcl
set_pine_mp_halo -halo_x 8 -halo_y 8 -min_channel {16 10}
set_pine_mp_halo -halo_x 8 -halo_y 8 -min_channel 0   ;# no channel rule
```

See [Minimum channel width](#minimum-channel-width) for what the rule does.

### set_pine_mp_row_alignment
Keep every hard macro's bottom edge on the standard-cell row grid (default: enabled).

```tcl
set_pine_mp_row_alignment -enable
set_pine_mp_row_alignment -disable
```

The grid is read from the rows: its origin is the lowest row and its pitch is the shortest site height. If the rows don't share one grid, PineMP warns (`PNE-0148`) and skips alignment. Macros whose height is a whole number of rows then have both horizontal edges on the grid; ariane's fakeram45_256x16 is 95 rows tall. See [Row alignment](#row-alignment) for why this matters.

### set_pine_mp_macro_halo
Set a per-macro halo override in microns. Two values are treated as symmetric horizontal/vertical padding; four values set left/bottom/right/top.

```tcl
set_pine_mp_macro_halo -macro_name MACRO1 -halo {10 10}
set_pine_mp_macro_halo -macro_name MACRO2 -halo {5 10 5 10}
```

### set_pine_mp_anchoring
Enable or disable four-corner anchoring (default: enabled).

```tcl
set_pine_mp_anchoring -enable
set_pine_mp_anchoring -disable
```

A B*-tree naturally compacts the macro cluster toward the bottom-left corner, which biases wirelength regardless of where the IO pins sit. After each SA run, PineMP re-evaluates the packed layout compacted toward each of the four core corners (bottom-left, bottom-right, top-left, top-right — obtained by reflecting the packing, following Chen & Chang, "Modern floorplanning based on B*-tree and fast simulated annealing") and keeps the corner with the lowest cost. This is cheap (no re-packing of the arrangement, only a coordinate reflection) and aligns the cluster with the fixed IO pins. Disable it to reproduce the classic bottom-left-only behavior. With [corner subtrees](#set_pine_mp_corner_subtrees), the four candidates are reflections of the whole layout: flipping in x swaps the left and right corner subtrees and mirrors the centered one.

### set_pine_mp_adaptive_weighting
Enable or disable adaptive IO/internal weight scheduling based on the current IO wirelength proportion.

```tcl
set_pine_mp_adaptive_weighting -enable
set_pine_mp_adaptive_weighting -disable
```

When enabled, PineMP computes the IO weight ratio from the current wirelength as
`(WLio * 100) / TWL` and uses that proportion to update the internal/IO weights dynamically during iteration.

### Shared halos

A halo is the area around a macro kept for buffers and the PDN, and `pine_mp` writes it out as a soft blockage. No other macro may be placed in it, but the halos of two neighbouring macros may overlap. Facing bodies therefore sit `max(h1, h2)` apart, the larger of their two facing halos, rather than `h1 + h2`. Each macro still keeps its full halo; the overlap belongs to both.

Node coordinates still give the halo box, with the body offset by the left and bottom halos, so node widths and heights, the cluster footprint and the outline penalty are unchanged. What changes:

1. **Packing.** A left child is placed against its neighbour's body rather than its halo box, so their facing halos overlap. The contour keeps two skylines, one of the halo boxes and one of the bodies. A new macro's body must clear the first and its halo the second. Over a strip covered only by the new macro's halo, the halo skyline keeps the higher of the old and new halo tops; otherwise a macro packed later could drop its body into a taller neighbour's halo.
2. **Cost.** The overlap penalty counts the area of each body inside another macro's halo box, both ways, instead of the overlap of the halo boxes. Halos still may not overlap pads, fixed instances or blockages.
3. **Slack.** [Slack moves](#set_pine_mp_slack_moves) measure how far a body can shift before it reaches a neighbour's halo, or brings its own halo onto the neighbour's body.
4. **Legalization.** The channel legalizer pushes along any macro whose body a move would bring into a halo, or whose halo it would bring onto a body. Its stop positions keep both halos clear.

With all halos zero, the packing is the same as before.

The [minimum channel](#minimum-channel-width) is a separate rule and still applies. It defaults to twice the halo, and two macros may abut only where neither facing side has a halo. So macros with the default halo still end up a full channel (two halos) apart, with their halos just touching. Halos actually overlap where the minimum channel is narrower than the two facing halos together: with per-macro halos larger than the default, or with a smaller `-min_channel`. Setting `-min_channel` to the halo lets default-halo macros share it fully (10 µm apart with 10 µm halos). But the ariane133/136 PDN failures came from channels exactly that narrow, so check `pdngen` and `analyze_power_grid` before relying on it.

### Minimum channel width

Pin-aware halos pad only the sides of a macro that carry signal pins. When a pin side faces a pinless side, or, since [halos are shared](#shared-halos), another pin side, the two macros end up exactly one halo apart, and B*-tree packing can also leave arbitrary gaps. Such slivers are filled with standard-cell rows that the rest of the flow cannot power. `tapcell` cuts rows 2 µm from every macro, and the PDN straps (metal4 at a 56 µm pitch on nangate45) may miss a narrow channel or stop just short of its last rail. The followpins then float, `pdngen` reports `PDN-0178 Remaining channel`, and `analyze_power_grid` fails with `PSM-0069`. Which macros land a few microns apart changes with every input, so the failure looked intermittent (ariane133/136: 8–10 µm channels with 8–10 µm halos).

The rule: after `pine_mp`, every gap between facing hard macros, and between a hard macro and the placement core edge, is either 0 (abutted) or at least the minimum channel for that axis. Two macros face each other when their projections overlap on the other axis. Side by side, they also face when only their corners meet: `tapcell` cuts every row within 2 µm of a macro, so if their vertical spans come within twice that, rounded out to whole rows (5.6 µm on nangate45), some row is cut from both ends and leaves a sliver between them. Gaps are measured between macro bodies and ignore halos. Soft macros are ignored because they are dropped after placement. Fixed macros are obstacles the rule applies to, but they are never moved, and gaps between two fixed macros are not checked. The rule is enforced at three levels:

1. **Packing.** The B*-tree contour keeps every vertical gap, between stacked macros and to the anchored horizontal core edge, either 0 or a full channel. A macro against the anchored vertical core edge is pushed off it when its halo would leave a sliver, and a macro placed beside its left neighbour in the tree is pushed out to a full channel when their shared halo would leave one. Corner anchoring is taken into account: the reflection moves the halo box but the body keeps its offset inside it, so which halo faces which neighbour depends on the corner.
2. **Annealing.** Packing only sees a macro's left neighbour in the tree, so horizontal gaps to other macros can still be narrow. `CostEvaluator` adds `0.05 × Σ (min − gap) / min` to the cost, so a one-halo sliver costs about 2.5% of the normalized wirelength. Flipping macros so that pinless sides face each other lets them abut and avoids the channel.
3. **Legalization.** Just before the final placement is written (`PNE-0050`), each remaining channel is either abutted (only if neither facing side has a halo) or widened to the minimum, by moving one side or both. If no channel can be fixed that way, the legalizer tries moving a side further, to its next abutment or full-channel position. As a last resort it slides a side sideways until the two no longer face each other. Moves are rigid pushes: any macro whose gap to a moved macro would shrink, whose halo the moved body would enter, or whose body would enter the moved halo, moves along with it. A move is kept only if no body enters another macro's halo, no halo overlaps an obstacle or leaves the core, and it strictly reduces the narrow channels (by count, then by total deficit). The result is deterministic. The minimum is rounded up to the manufacturing grid, so all shifts are multiples of it.

`PNE-0142` logs the minimum channel in effect, and `PNE-0140` logs `Narrow channels: N found, F fixed, R remaining`. Every channel that could not be fixed gets a `PNE-0141` warning naming both macros (or the core edge) and the channel coordinates. With `-uniform` halos, facing macros are always at least one halo apart, so a minimum channel no wider than the halo only affects the core edges. `-min_channel 0` together with `set_pine_mp_row_alignment -disable` turns both rules off. The placements still differ from those before either rule existed, because halos are now [shared](#shared-halos).

### Row alignment

The channel rule alone did not stop the PDN failures. After it, ariane133 (10 µm halos, 8 partitions) still failed `PSM-0069` on VDD. The channels were now about 22 µm wide, and `pdngen` repaired them with a metal4 strap down the middle (`PDN-0178` was gone). The rail that failed was the one nearest the macro that closed the channel off: the 18 µm VDD segment at y = 1085.0 µm, the bottom rail of the lowest row in the channel. `pdngen` had built its only via and then removed it (`PDN-0195 Removing 9 via(s) between metal1 and metal4 at (1194.855, 1085.075)`).

The cause is that two tools keep their distance from a macro edge in different ways:

- **The rows** start where `tapcell` puts them: 2 µm from the macro, rounded *outward to the next row boundary* (rows exist only on the 1.4 µm row grid).
- **The metal4 straps** end where `pdngen` cuts them: at the macro's PDN halo (2 µm) plus the strap's spacing (0.27 µm for a 0.48 µm strap on nangate45), measured *exactly* from the macro edge.

So the gap between the first rail and the strap end depends on where the macro edge sits relative to the row grid. Call that offset *r*, how far the edge is above the row boundary below it. For a macro top edge with rows above it, the first rail is at the next row boundary at least 2 µm above the edge. It is 0.17 µm wide, and the strap starts 2.27 µm above the edge:

| r (µm) | first rail above the edge | strap start | result |
|---|---|---|---|
| 0, 0.2, 0.4, 1.0, 1.2 | 2.4–3.2 µm | 2.27 µm | rail fully on the strap |
| 0.6 | 2.2 µm (2.115–2.285) | 2.27 µm | 0.015 µm overlap: the via is built, then removed (`PDN-0195`) |
| 0.8 | 2.0 µm (1.915–2.085) | 2.27 µm | no overlap, no via |

Below a macro's bottom edge the rows run downward, and the same two offsets, 0.6 and 0.8, fail. A long rail survives, because other straps further along still connect it. A rail in a channel closed off by a macro has only the one strap, so it floats. PineMP's halos and channels were not whole rows (10 µm = 7 rows + 0.2, 20 µm = 14 rows + 0.4), so macro edges landed at offsets that are multiples of 0.2 µm, two of which fail. The ariane133 log matches r = 0.6: the removed via is centred at 1085.075 µm, inside the 1085.07–1085.085 overlap, and a removed VSS via at y = 498.475 µm sits the same distance above its rail. Whether a closing macro hit 0.6 or 0.8 changed with every arrangement, which is why the failure came and went with the partition count.

With the edge on the grid (r = 0), the first rail is 2.8 µm from the macro on both sides and the strap covers it with 0.45 µm to spare. PineMP enforces this at two levels:

1. **Packing.** `pack()` raises each hard macro to the next position where its body, once the cluster is reflected to its anchor corner, starts on a row. Unlike the channel spacing, this stays on in the fallback packing, so the SA never sees unaligned states. It costs at most one row per macro in the stack.
2. **Legalization.** The same pass that fixes channels counts off-grid macros as violations and moves them to the row below or above. It never accepts a channel fix that takes a macro off its row, and it widens y channels in whole rows. `PNE-0145` logs the grid, `PNE-0146` logs `Macros off the row grid: N found, F fixed, R remaining`, and every macro left off the grid gets a `PNE-0147` warning.

The safe offsets in the table assume the ORFS nangate45 setup: a 2 µm `tapcell` row halo and a 2 µm macro halo on the PDN grid. With a different setup, check that the first rail beyond a row-aligned macro (the row halo rounded up to whole rows) sits clear of where the straps end (the PDN halo plus the strap spacing).

### Corner neighbours

Row alignment fixed the via removals, but ariane133 (8 partitions) then failed `PSM-0069` on VSS in a pocket at x 1098–1124 µm, y 1103–1143 µm, enclosed by about five macros. pdngen reported it as `PDN-0178 Remaining channel ... on metal1 for nets: VSS`. VDD survived only because a regular VDD strap happened to fall inside. The endcap cells show the pocket was a staircase:

- At the top-left, rows 810–812 were a 5.7 µm sliver between one macro ending at y = 1142.4 µm and another starting at y = 1142.4 µm, 10 µm to its right.
- At the bottom-right, rows 783–784 formed a notch between two macros whose corners met at y = 1101.8 µm.

Neither pair counted as facing, because their spans overlapped by zero. Yet `tapcell` cut those rows from both ends. pdngen's channel repair places one straight strap that must clear obstructions over the full height of the channel's rails. Here the bottom-left macro blocked x < 1108.4 µm near the bottom and the top-right macro blocked x > 1104.3 µm near the top, so no position worked. The sliver and notch rails could only get short stub straps, which reach no metal7. With corner neighbours counted as facing, the 10 µm pair is widened to a full channel. That leaves a strap position clear of both macros, and the sliver becomes a normal row segment.

`PNE-0142` reports the vertical distance within which side-by-side macros count as facing. It assumes `tapcell`'s default 2 µm row halo.

## Example Usage

```tcl
# Load design
read_lef tech.lef
read_def design.def

# Configure PineMP
set_pine_mp_iterations 5
set_pine_mp_initial_weights -internal_weight 0.8 -io_weight 0.2
set_pine_mp_final_weights -internal_weight 0.5 -io_weight 0.5
set_pine_mp_pin_strategy connectivity
set_pine_mp_halo -halo_x 10 -halo_y 10 -uniform
set_pine_mp_adaptive_weighting -enable

# Run placement
pine_mp

# Save result
write_def placed_design.def
```

## Soft Macros

PineMP supports soft macro regions built from netlist partitions. Use `set_pine_mp_soft_macros` to enable co-placement of soft macro clusters alongside hard macros.

```tcl
set_pine_mp_soft_macros -utilization 0.7 -aspect_ratio 1.0
```

Soft macros are represented virtually in the B*-tree and are not converted to DB blockages; their final positions are reported separately.

### Standard Cell Seeding

```tcl
set_pine_mp_stdcell_seeding -enable
set_pine_mp_stdcell_seeding -disable
```

Enabled by default, and only active when soft macros are enabled. Every
standard cell is moved to the center of the soft macro its partition produced
and marked `PLACED`. All the cells of a partition land on the same point: this
is not a legal placement, it is a *seed* for pin placement and global
placement.

The point of the seed is that `ppl` and `gpl` discard the partitioning
otherwise. Both read an instance's DB location only when the instance is
placed, and put every other instance at the center of the die (`ppl`) or of
the core (`gpl`, whose `GPL-0051` reports the split). Without seeding, the
soft macro arrangement PineMP optimizes is lost for the standard cells and
survives only through the hard macro positions. `mpl` does the same thing
with its hierarchy clusters (`HierRTLMP::setTemporaryStdCellLocation`).

PineMP seeds the cells before every `ppl` call, so the IO pins follow the soft
macros of the current iteration, and once more at the end, which `PNE-0094`
reports. The seed is the soft macro center, the same point `CostEvaluator`
uses for every member pin, so the wirelength PineMP reports is the wirelength
of the configuration handed to `gpl`. The cost model maps member cells to
their soft macro before it looks at their placement, so seeding changes what
`ppl` sees but never the cost.

After each iteration's `ppl` call, `PNE-0066` reports how the IO wirelength
changed, and under `-style spread` `PNE-0067` reports the soft IO wirelength
and the hard extension. With the seed, `ppl` and the cost model see nearly the
same pin positions, so the IO wirelength should usually drop. Without it,
`ppl` aims every pin that reaches the standard cells at the die center.

Instances that are fixed (`LOCKED`, `FIRM`, `COVER`) are left untouched.
Disable seeding to hand `ppl` and `gpl` an unplaced netlist, as before
seeding existed.

### Internal Recursive Partitioning

By default `pine_mp` partitions the standard cells itself using recursive bisection (divide and conquer): the whole netlist is split in two with the par module's min-cut engine, then the largest partition is repeatedly split in half until the requested number of partitions and/or maximum partition size is reached. This is faster than a single flat k-way `triton_part` call, avoids the external call in the flow, and preserves the bisection hierarchy for future hierarchical placement.

```tcl
set_pine_mp_partitioning [-num_partitions num] [-max_area_fraction fraction] \
                         [-min_cells num] [-seed seed] [-external]
```

- `-num_partitions` — split until at least this many partitions exist (count target).
- `-max_area_fraction` — keep splitting any partition whose soft-macro footprint (cell area / utilization) would exceed this fraction of the core area. Useful for designs like ariane133/136 where soft macros would otherwise dwarf the hard macros. Takes precedence over the count target.
- `-min_cells` — never split groups below this cell count (default 50).
- `-seed` — partitioner seed (default 1).
- `-external` — skip internal partitioning and consume pre-existing `partition_id` properties (e.g. from a manual `triton_part` run).

When neither `-num_partitions` nor `-max_area_fraction` is given, 10 partitions are created (matching the previous external flow default).

```tcl
# Example: cap each soft macro at 10% of the core, at least 16 partitions
set_pine_mp_partitioning -num_partitions 16 -max_area_fraction 0.1
set_pine_mp_soft_macros -utilization 0.7
pine_mp
pine_mp_report_partition_tree
```

`pine_mp_report_partition_tree` prints the bisection hierarchy of the last run (internal nodes and leaf partitions with cell counts and areas).

## Architecture

The implementation is modular and research-friendly:

- **BStarTree** - B*-Tree representation for macro placement: one tree, or a forest of corner subtrees around a centered soft one
- **ChannelLegalizer** - Narrow-channel and off-row detection, SA penalty and final legalization
- **CostEvaluator** - Weighted HPWL cost computation with net classification
- **WeightScheduler** - Adaptive weight management for iterative optimization
- **PinAssigner** - Pin assignment engine with multiple strategies
- **SimulatedAnnealing** - SA optimizer with configurable perturbations
- **PineMP** - Main coordinator orchestrating the iterative loop

## Research Notes

This is experimental research code designed for flexibility and easy experimentation with different weight schedules, pin assignment strategies, and cost functions.
