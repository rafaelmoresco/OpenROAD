# PineMP Basic Test
# Tests the basic iterative co-optimization flow

source "helpers.tcl"
read_lef Nangate45/Nangate45.lef
read_def one_macro.def

# Configure PineMP with basic settings
set_pine_mp_iterations 3
set_pine_mp_initial_weights -internal_weight 0.8 -io_weight 0.2
set_pine_mp_final_weights -internal_weight 0.5 -io_weight 0.5
set_pine_mp_sa_params -initial_temp 500.0 -cooling_rate 0.9 -max_iterations 1000
set_pine_mp_pin_strategy connectivity

# Run PineMP
pine_mp

set def_file [make_result_file pine_basic.def]
write_def $def_file
diff_file pine_basic.defok $def_file
