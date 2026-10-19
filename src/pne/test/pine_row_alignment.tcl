# Hard macro edges must sit on the std-cell row grid.  The PDN straps stop a
# fixed distance from each macro while the rows next to it stop on the row
# grid, so an off-grid edge can leave that row's rail only partly under the
# strap and unconnected.
source "pine_test_utils.tcl"

pne_load_design "narrow_channel.def" "narrow_channel.lef"

# Only the row grid, no channel rule.  A half-row vertical halo puts the
# macro against the anchored core edge off the 1.4 um rows, whichever corner
# is chosen (see pine_row_alignment_off).
set_pine_mp_halo -halo_x 5 -halo_y 0.7 -uniform -min_channel 0
set_pine_mp_iterations 0

pne_run_and_save pine_row_alignment
pne_assert_macros_inside_core

set off_grid [pne_off_grid_macros]
if { [llength $off_grid] != 0 } {
  error "macros off the row grid: $off_grid"
}

puts "pass"
