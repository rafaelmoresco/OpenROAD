# With corner subtrees, MACRO_1 is the root of its corner subtree, so it ends
# flush with a core corner (a top corner lowers it by up to one row to stay
# on the row grid), and the two equal soft macros form the centered subtree,
# so the midpoint of their centers, where pine_mp seeds the cells, is the
# core center.
source "pine_test_utils.tcl"

pne_load_design "soft_clusters.def" "macro_only.lef"
pne_set_quick_defaults
set_pine_mp_soft_macros
set_pine_mp_partitioning -external
set_pine_mp_objective -style spread
set_pine_mp_corner_subtrees -enable -soft_overlap_weight 0.2

pne_run_and_save pine_corner_subtrees
pne_assert_macros_inside_core

set block [ord::get_db_block]
set core [$block getCoreArea]
set row_pitch 2800

set macro [$block findInst MACRO_1]
if { [$macro getPlacementStatus] != "LOCKED" } {
  error "MACRO_1 is [$macro getPlacementStatus] after pine_mp, expected LOCKED"
}
set box [$macro getBBox]
set dx [expr { min([$box xMin] - [$core xMin], [$core xMax] - [$box xMax]) }]
set dy [expr { min([$box yMin] - [$core yMin], [$core yMax] - [$box yMax]) }]
if { $dx != 0 || $dy < 0 || $dy >= $row_pitch } {
  error "MACRO_1 at ([$box xMin], [$box yMin]) is $dx, $dy from a core corner"
}

proc center { inst_name } {
  set box [[[ord::get_db_block] findInst $inst_name] getBBox]
  return [list [expr { ([$box xMin] + [$box xMax]) / 2 }] \
            [expr { ([$box yMin] + [$box yMax]) / 2 }]]
}
lassign [center c0_0] x0 y0
lassign [center c1_0] x1 y1
set cx [expr { ([$core xMin] + [$core xMax]) / 2 }]
set cy [expr { ([$core yMin] + [$core yMax]) / 2 }]
set off [expr { abs(($x0 + $x1) / 2 - $cx) + abs(($y0 + $y1) / 2 - $cy) }]
if { $off > 2000 } {
  error "soft macro centers ($x0, $y0), ($x1, $y1): $off from the core center"
}

puts "pass"
