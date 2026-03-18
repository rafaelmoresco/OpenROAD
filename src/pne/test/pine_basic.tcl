# PineMP basic flow with explicit settings.
source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

set_pine_mp_iterations 3
set_pine_mp_initial_weights -internal_weight 0.8 -io_weight 0.2
set_pine_mp_final_weights -internal_weight 0.5 -io_weight 0.5
set_pine_mp_sa_params -initial_temp 300.0 -cooling_rate 0.9 -max_iterations 250
set_pine_mp_pin_strategy connectivity

pne_run_and_save pine_basic

puts "pass"
