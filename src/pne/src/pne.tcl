# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021-2025, The OpenROAD Authors

sta::define_cmd_args "pine_mp" { [-num_threads num_threads] }

proc pine_mp { args } {
  sta::parse_key_args "pine_mp" args \
    keys {-num_threads} \
    flags {}

  sta::check_argc_eq0 "pine_mp" $args

  #
  # Check for valid design
  if { [ord::get_db_block] == "NULL" } {
    utl::error PNE 3 "No block found for PineMP placement."
  }

  # Set default parameters
  set num_threads 1

  if { [info exists keys(-num_threads)] } {
    set num_threads $keys(-num_threads)
  }

  if {
    ![pne::pine_mp_cmd $num_threads]
  } {
    return false
  }

  return true
}

sta::define_cmd_args "pine_mp_debug" { [-debug] }

proc pine_mp_debug { args } {
  sta::parse_key_args "pine_mp_debug" args \
    keys {} \
    flags {-debug}

  set debug [info exists flags(-debug)]

  pne::set_debug_cmd $debug
}

sta::define_cmd_args "set_pine_mp_iterations" { num_iterations }

proc set_pine_mp_iterations { num_iterations } {
  pne::set_num_iterations_cmd $num_iterations
}

sta::define_cmd_args "set_pine_mp_initial_weights" \
  { -internal_weight internal_weight -io_weight io_weight }

proc set_pine_mp_initial_weights { args } {
  sta::parse_key_args "set_pine_mp_initial_weights" args \
    keys {-internal_weight -io_weight} \
    flags {}

  if { ![info exists keys(-internal_weight)] } {
    utl::error PNE 100 "Missing required argument -internal_weight"
  }
  
  if { ![info exists keys(-io_weight)] } {
    utl::error PNE 101 "Missing required argument -io_weight"
  }

  pne::set_initial_weights_cmd $keys(-internal_weight) $keys(-io_weight)
}

sta::define_cmd_args "set_pine_mp_final_weights" \
  { -internal_weight internal_weight -io_weight io_weight }

proc set_pine_mp_final_weights { args } {
  sta::parse_key_args "set_pine_mp_final_weights" args \
    keys {-internal_weight -io_weight} \
    flags {}

  if { ![info exists keys(-internal_weight)] } {
    utl::error PNE 102 "Missing required argument -internal_weight"
  }
  
  if { ![info exists keys(-io_weight)] } {
    utl::error PNE 103 "Missing required argument -io_weight"
  }

  pne::set_final_weights_cmd $keys(-internal_weight) $keys(-io_weight)
}

sta::define_cmd_args "set_pine_mp_sa_params" \
  { [-initial_temp temp] [-cooling_rate rate] [-max_iterations iterations] }

proc set_pine_mp_sa_params { args } {
  sta::parse_key_args "set_pine_mp_sa_params" args \
    keys {-initial_temp -cooling_rate -max_iterations} \
    flags {}

  # Defaults
  set initial_temp 1000.0
  set cooling_rate 0.95
  set max_iterations 10000

  if { [info exists keys(-initial_temp)] } {
    set initial_temp $keys(-initial_temp)
  }
  
  if { [info exists keys(-cooling_rate)] } {
    set cooling_rate $keys(-cooling_rate)
  }
  
  if { [info exists keys(-max_iterations)] } {
    set max_iterations $keys(-max_iterations)
  }

  pne::set_sa_params_cmd $initial_temp $cooling_rate $max_iterations
}

sta::define_cmd_args "set_pine_mp_pin_strategy" { strategy }

proc set_pine_mp_pin_strategy { strategy } {
  # Valid strategies: uniform, connectivity, random, hungarian
  pne::set_pin_strategy_cmd $strategy
}

namespace eval pne {
  # Internal utility functions can be added here
}
