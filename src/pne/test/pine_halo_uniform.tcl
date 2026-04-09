# Test uniform halo (all four sides, ignoring pin locations).
source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

# Use -uniform to apply halo on all sides regardless of pin locations.
set_pine_mp_halo -halo_x 5 -halo_y 5 -uniform

pne_set_quick_defaults
pne_run_and_save pine_halo_uniform

pne_assert_macros_inside_core

puts "pass"
