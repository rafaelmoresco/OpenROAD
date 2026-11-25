# The std cells are seeded at their soft macro centers before every ppl call,
# so io_0 and io_1, which reach one cluster each, end up about as close to it
# as the die edge allows.  Unseeded, ppl aims both at the die center.
source "pine_test_utils.tcl"

pne_load_design "soft_clusters.def" "macro_only.lef"
pne_set_quick_defaults
set_pine_mp_soft_macros
set_pine_mp_partitioning -external

pne_run_and_save pine_ppl_seeding
pne_assert_macros_inside_core

set block [ord::get_db_block]
set die [$block getDieArea]
# 20 um: ppl's sections and the pin offset inside a cell keep a pin off the
# exact point of the die edge nearest its cluster.
set tolerance 40000

foreach cluster { 0 1 } {
  set cell [[$block findInst "c${cluster}_0"] getBBox]
  set cx [expr { ([$cell xMin] + [$cell xMax]) / 2 }]
  set cy [expr { ([$cell yMin] + [$cell yMax]) / 2 }]
  set edge_dist [expr { min($cx - [$die xMin], [$die xMax] - $cx,
                            $cy - [$die yMin], [$die yMax] - $cy) }]

  set pin [[$block findBTerm "io_$cluster"] getBBox]
  set px [expr { ([$pin xMin] + [$pin xMax]) / 2 }]
  set py [expr { ([$pin yMin] + [$pin yMax]) / 2 }]
  set dist [expr { abs($px - $cx) + abs($py - $cy) }]

  if { $dist > $edge_dist + $tolerance } {
    error "io_$cluster at ($px, $py) is $dist from its cluster at\
      ($cx, $cy), which is $edge_dist from the die edge"
  }
}

puts "pass"
