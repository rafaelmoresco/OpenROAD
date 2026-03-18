source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"
set_pine_mp_iterations 1
set_pine_mp_initial_weights -internal_weight 0.8 -io_weight 0.2
set_pine_mp_final_weights -internal_weight 0.6 -io_weight 0.4
set_pine_mp_sa_params -initial_temp 200.0 -cooling_rate 0.9 -max_iterations 100

foreach strategy {uniform connectivity random hungarian unknown_strategy} {
  set_pine_mp_pin_strategy $strategy
  if { ![pine_mp -num_threads 1] } {
    error "pine_mp failed for strategy $strategy"
  }
}

puts "pass"
