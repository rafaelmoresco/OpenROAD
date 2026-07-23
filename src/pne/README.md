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
   - Cost function uses weighted HPWL: `α * Internal_WL + β * IO_WL + penalties`
   - Initially α (internal weight) is high, β (IO weight) is low
   - After SA convergence: reassign pins based on current placement
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

Following Adya & Markov (ParquetFP, "Fixed-outline floorplanning: enabling hierarchical design"), each block is assigned a *slack* in x and y — how far it can shift toward the right / top before hitting a neighbour or the used extent. Zero-slack blocks are *critical*: they set the floorplan's width / height. A fraction `p` of SA moves (default 0.5) are biased to pick a critical block in whichever dimension is currently most over the core outline and relocate it toward a block with spare room (high slack); the rest stay uniform-random for ergodicity. This packs the violating dimension faster than blind moves and reduces out-of-bounds macros. Slack is recomputed once per temperature step. Set `-probability 0.0` (or `-disable`) to recover pure uniform-random moves.

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

### set_pine_mp_macro_halo
Set a per-macro halo override in DBU units. Two values are treated as symmetric horizontal/vertical padding; four values set left/bottom/right/top.

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

A B*-tree naturally compacts the macro cluster toward the bottom-left corner, which biases wirelength regardless of where the IO pins sit. After each SA run, PineMP re-evaluates the packed layout compacted toward each of the four core corners (bottom-left, bottom-right, top-left, top-right — obtained by reflecting the packing, following Chen & Chang, "Modern floorplanning based on B*-tree and fast simulated annealing") and keeps the corner with the lowest cost. This is cheap (no re-packing of the arrangement, only a coordinate reflection) and aligns the cluster with the fixed IO pins. Disable it to reproduce the classic bottom-left-only behavior.

### set_pine_mp_adaptive_weighting
Enable or disable adaptive IO/internal weight scheduling based on the current IO wirelength proportion.

```tcl
set_pine_mp_adaptive_weighting -enable
set_pine_mp_adaptive_weighting -disable
```

When enabled, PineMP computes the IO weight ratio from the current wirelength as
`(WLio * 100) / TWL` and uses that proportion to update the internal/IO weights dynamically during iteration.

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

- **BStarTree** - B*-Tree representation for macro placement
- **CostEvaluator** - Weighted HPWL cost computation with net classification
- **WeightScheduler** - Adaptive weight management for iterative optimization
- **PinAssigner** - Pin assignment engine with multiple strategies
- **SimulatedAnnealing** - SA optimizer with configurable perturbations
- **PineMP** - Main coordinator orchestrating the iterative loop

## Research Notes

This is experimental research code designed for flexibility and easy experimentation with different weight schedules, pin assignment strategies, and cost functions.
