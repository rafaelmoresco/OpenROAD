# Neighbouring halos may overlap, but no macro may sit inside another macro's
# halo.  With the channel rule and row alignment off, the initial packing puts
# every neighbour one shared halo away.
source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"
set_pine_mp_halo -halo_x 5 -halo_y 5 -uniform -min_channel 0
set_pine_mp_row_alignment -disable
set_pine_mp_iterations 0

pne_run_and_save pine_halo_shared
pne_assert_macros_inside_core

lassign [pne_halo_sharing] intrusions shared
if { [llength $intrusions] != 0 } {
  error "macros inside another macro's halo: $intrusions"
}
if { $shared == 0 } {
  error "expected neighbouring halos to overlap"
}

puts "pass"
