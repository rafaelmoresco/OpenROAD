# pine_min_channel with the rule (and row alignment) disabled: the one-halo
# slivers stay, which shows the scenario does produce narrow channels and
# that 0 turns it off.
source "pine_test_utils.tcl"

pne_load_design "narrow_channel.def" "narrow_channel.lef"
set_pine_mp_halo -halo_x 5 -halo_y 5 -min_channel 0
set_pine_mp_row_alignment -disable
set_pine_mp_iterations 0

pne_run_and_save pine_min_channel_off
pne_assert_macros_inside_core

set min_channel [ord::microns_to_dbu 10]
if { [llength [pne_narrow_channels $min_channel $min_channel]] == 0 } {
  error "expected the unconstrained packing to leave narrow channels"
}

puts "pass"
