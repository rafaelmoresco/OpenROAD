# PineMP Current Implementation
## B*-Tree and Simulated Annealing Cost Function

---

## Slide 1 — PineMP Submodule: Current Contents

- **Purpose**: Pin-aware iterative macro placement with co-optimization of placement and pin assignment.
- **Core classes**:
  - `PineMP`: main coordinator for iterative flow
  - `BStarTree`: floorplan representation and packing
  - `SimulatedAnnealing`: perturb/accept optimization loop
  - `CostEvaluator`: weighted HPWL + penalties
  - `WeightScheduler`: iteration-wise weight updates
  - `PinAssigner`: uniform/connectivity/random strategies
- **User-facing Tcl commands**:
  - `pine_mp`
  - `set_pine_mp_iterations`
  - `set_pine_mp_initial_weights`, `set_pine_mp_final_weights`
  - `set_pine_mp_sa_params`
  - `set_pine_mp_pin_strategy`

**Speaker notes:**
- Emphasize modularity: each optimization concern is isolated in its own class.
- Mention this is an iterative co-optimization loop, not a one-shot placement.
- Optional caveat: Hungarian pin strategy is declared but currently falls back to connectivity.

---

## Slide 2 — B*-Tree: Current Implementation

- **Data structure**:
  - Each macro is a `BStarNode` with `(parent, left, right)`, orientation, and packed `(x, y)`.
- **Initial tree build**:
  - Macros are inserted in a simple right-child chain from the root.
- **Packing algorithm (contour-based)**:
  - Root starts at origin.
  - Left child is placed to the right of parent (`x + width`).
  - Right child is placed at same `x` and compacted upward by contour height.
  - Bounding box (`width`, `height`, `area`) is recomputed after packing.
- **SA perturbation operators on tree**:
  - `swapNodes(id1, id2)`
  - `rotateNode(id)` (R0→R90→R180→R270)
  - `moveNode(id, new_parent_id, as_left_child)`
- **State management**:
  - `save()` / `restore()` snapshot structure, orientation, and coordinates.
  - `applyPlacement()` writes macro locations/status to OpenDB.

**Speaker notes:**
- Highlight that left/right semantics match common B*-Tree packing intuition (rightward and upward compaction).
- Mention this implementation favors simplicity for experimentation.
- Optional caveat: move/remove reinsertion logic is simplified compared to production-grade legalizers.

---

## Slide 3 — SA Cost Function: Current Implementation

- **Objective in SA**:

  \[
  C = C_{wl} + w_{ov} \cdot C_{overlap} + w_{out} \cdot C_{outline}
  \]

  with

  \[
  C_{wl} = w_{int} \cdot WL_{internal} + w_{io} \cdot WL_{io}
  \]

- **Net classification for wirelength**:
  - `INTERNAL`: macro-only nets
  - `IO`: nets touching block terminals (BTerms)
  - `MIXED`: split 50/50 between internal and IO contributions
  - Power/ground/clock nets are skipped.
- **Penalty terms**:
  - `C_overlap`: pairwise macro overlap area
  - `C_outline`: overflow penalty when packed width/height exceed die limits
  - Default penalty weights are high (`w_{ov}=1e6`, `w_{out}=1e6`).
- **Annealing acceptance**:
  - Always accept improvements (`\Delta C < 0`)
  - Otherwise accept with probability `exp(-\Delta C / T)`.
- **Weight schedule across outer iterations**:
  - Default linear transition from `(w_int, w_io) = (0.8, 0.2)` to `(0.5, 0.5)`.

**Speaker notes:**
- Explain that high overlap/outline weights push feasibility strongly.
- Clarify that internal-vs-IO balance is intentionally shifted over iterations.
- Optional caveat: current implementation tracks best cost but does not maintain a separate best-geometry snapshot.
