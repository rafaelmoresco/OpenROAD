# pine_row_alignment with the row grid disabled: the half-row halo leaves the
# macro against the anchored core edge off the rows, which shows the scenario
# needs the alignment.
source "pine_test_utils.tcl"

pne_load_design "narrow_channel.def" "narrow_channel.lef"
set_pine_mp_halo -halo_x 5 -halo_y 0.7 -uniform -min_channel 0
set_pine_mp_row_alignment -disable
set_pine_mp_iterations 0

pne_run_and_save pine_row_alignment_off
pne_assert_macros_inside_core

if { [llength [pne_off_grid_macros]] == 0 } {
  error "expected the unaligned packing to leave macros off the row grid"
}

puts "pass"
