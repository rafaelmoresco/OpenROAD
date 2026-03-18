source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

# Exercise default arguments for SA and debug command plumbing.
set_pine_mp_iterations 3
set_pine_mp_initial_weights -internal_weight 0.7 -io_weight 0.3
set_pine_mp_final_weights -internal_weight 0.5 -io_weight 0.5
set_pine_mp_sa_params
pine_mp_debug -debug

set def_file [pne_run_and_save pine_api_defaults]
pne_assert_def_origins_within_die $def_file

puts "pass"
