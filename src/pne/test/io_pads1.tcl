source "pine_test_utils.tcl"

set def_file [pne_run_scenario "io_pads1" "io_pads1.def" "macro_only.lef" "dummy_pads.lef"]
pne_assert_macros_clear_io_pads
pne_assert_def_origins_within_die $def_file

puts "pass"
