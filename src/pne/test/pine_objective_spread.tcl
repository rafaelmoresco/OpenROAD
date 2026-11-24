# With -style spread, hard macros are pushed away from their IO pins.  The
# only macro in fixed_ios1 drives a fixed pin on the left edge, so it must
# end up right of the core center (min_wl puts it left of it).
source "pine_test_utils.tcl"

pne_load_design "fixed_ios1.def" "macro_only.lef"
pne_set_quick_defaults
set_pine_mp_objective -style spread -spread_weight 2.0

pne_run_and_save pine_objective_spread
pne_assert_macros_inside_core

set block [ord::get_db_block]
set core [$block getCoreArea]
set bbox [[$block findInst MACRO_1] getBBox]
set macro_x [expr { ([$bbox xMin] + [$bbox xMax]) / 2 }]
set core_x [expr { ([$core xMin] + [$core xMax]) / 2 }]
if { $macro_x <= $core_x } {
  error "MACRO_1 center x $macro_x is not right of the core center $core_x"
}

puts "pass"
