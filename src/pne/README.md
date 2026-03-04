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

### set_pine_mp_pin_strategy
Set the pin assignment strategy (uniform, connectivity, random, hungarian).

```tcl
set_pine_mp_pin_strategy connectivity
```

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

# Run placement
pine_mp

# Save result
write_def placed_design.def
```

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
