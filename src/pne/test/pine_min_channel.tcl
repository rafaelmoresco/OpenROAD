# Pin-aware halos pad only pin sides, so a pin side facing a pinless side
# leaves a one-halo sliver that std-cell rows and the PDN cannot use.  Every
# gap must end up abutted or at least the minimum channel wide.
source "pine_test_utils.tcl"

pne_load_design "narrow_channel.def" "narrow_channel.lef"

# Default minimum channel: 2x the halo, i.e. 10 um (11.2 um, whole rows, along
# y with the default row alignment).
set_pine_mp_halo -halo_x 5 -halo_y 5

# No annealing: the initial packing alone leaves narrow channels (see
# pine_min_channel_off), so this checks pack() and the final legalization.
set_pine_mp_iterations 0

pne_run_and_save pine_min_channel
pne_assert_macros_inside_core

set min_channel [ord::microns_to_dbu 10]
set row_cut [ord::microns_to_dbu 2.8]
set narrow [pne_narrow_channels $min_channel $min_channel $row_cut]
if { [llength $narrow] != 0 } {
  error "narrow channels left: $narrow"
}

set off_grid [pne_off_grid_macros]
if { [llength $off_grid] != 0 } {
  error "macros off the row grid: $off_grid"
}

puts "pass"
