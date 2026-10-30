# Full flow (multi-start annealing, anchoring, legalization) with separate
# x and y minimum channels.
source "pine_test_utils.tcl"

pne_load_design "narrow_channel.def" "narrow_channel.lef"
set_pine_mp_halo -halo_x 5 -halo_y 5 -min_channel {12 10}
pne_set_quick_defaults

pne_run_and_save pine_min_channel_sa 2
pne_assert_macros_inside_core

set narrow [pne_narrow_channels [ord::microns_to_dbu 12] [ord::microns_to_dbu 10] \
              [ord::microns_to_dbu 2.8]]
if { [llength $narrow] != 0 } {
  error "narrow channels left: $narrow"
}

set off_grid [pne_off_grid_macros]
if { [llength $off_grid] != 0 } {
  error "macros off the row grid: $off_grid"
}

puts "pass"
