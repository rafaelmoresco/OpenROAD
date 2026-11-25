# With hiding, ppl reads MACRO_1 at the die center in every call, so io_m,
# which reaches only MACRO_1, goes to the middle of the bottom or top die
# edge instead of next to MACRO_1's pin.  MACRO_1 must end up LOCKED again.
source "pine_test_utils.tcl"

pne_load_design "soft_clusters.def" "macro_only.lef"
pne_set_quick_defaults
set_pine_mp_ppl_macro_hiding -enable

pne_run_and_save pine_ppl_macro_hiding
pne_assert_macros_inside_core

set block [ord::get_db_block]
set status [[$block findInst MACRO_1] getPlacementStatus]
if { $status != "LOCKED" } {
  error "MACRO_1 is $status after pine_mp, expected LOCKED"
}

# The die center is 125000 from the bottom and top edges and 150000 from the
# sides.  MACRO_1's pin sits halfway up the macro, nearer a side.
set die [$block getDieArea]
set cx [expr { ([$die xMin] + [$die xMax]) / 2 }]
set cy [expr { ([$die yMin] + [$die yMax]) / 2 }]
set limit [expr { ([$die dx] + [$die dy]) / 4 }]

set pin [[$block findBTerm io_m] getBBox]
set px [expr { ([$pin xMin] + [$pin xMax]) / 2 }]
set py [expr { ([$pin yMin] + [$pin yMax]) / 2 }]
set dist [expr { abs($px - $cx) + abs($py - $cy) }]
if { $dist > $limit } {
  error "io_m at ($px, $py) is $dist from the die center, more than $limit"
}

puts "pass"
