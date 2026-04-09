# Test basic halo functionality with pin-aware placement.
source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

# Set a halo of 5 microns around pin sides.
set_pine_mp_halo -halo_x 5 -halo_y 5

pne_set_quick_defaults
pne_run_and_save pine_halo_basic

# Verify macros do not exceed the core area even with halos.
pne_assert_macros_inside_core

puts "pass"
