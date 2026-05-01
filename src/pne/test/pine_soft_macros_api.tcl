# Test: set_pine_mp_soft_macros and pine_mp_report_soft_macros API.
#
# This test verifies that the soft macro Tcl commands are properly wired up
# and that calling them on a design without prior partition_id properties
# emits the expected warning (no crash).

source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

# Configure soft macros with custom utilization and aspect ratio.
set_pine_mp_soft_macros -utilization 0.65 -aspect_ratio 1.2

# Run placement.  The design has no partition_id properties yet, so PineMP
# should warn and continue placing only the hard macros.
set_pine_mp_iterations 2
pne_set_quick_defaults

pine_mp

# Reporting soft macros before they are populated should not crash.
pine_mp_report_soft_macros

puts "pass"
