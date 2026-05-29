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
    utl::error PNE 99 "No block found for PineMP placement."
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

sta::define_cmd_args "set_pine_mp_adaptive_weighting" \
  { [-enable] [-disable] }

proc set_pine_mp_adaptive_weighting { args } {
  sta::parse_key_args "set_pine_mp_adaptive_weighting" args \
    keys {} \
    flags {-enable -disable}

  set adaptive 1
  if { [info exists flags(-disable)] } {
    set adaptive 0
  }

  pne::set_adaptive_weighting_cmd $adaptive
}

sta::define_cmd_args "set_pine_mp_halo" \
  { -halo_x halo_x -halo_y halo_y [-pin_aware] [-uniform] }

proc set_pine_mp_halo { args } {
  sta::parse_key_args "set_pine_mp_halo" args \
    keys {-halo_x -halo_y} \
    flags {-pin_aware -uniform}

  if { ![info exists keys(-halo_x)] } {
    utl::error PNE 110 "Missing required argument -halo_x"
  }

  if { ![info exists keys(-halo_y)] } {
    utl::error PNE 111 "Missing required argument -halo_y"
  }

  set halo_x [ord::microns_to_dbu $keys(-halo_x)]
  set halo_y [ord::microns_to_dbu $keys(-halo_y)]

  # Default is pin_aware unless -uniform is specified
  set pin_aware 1
  if { [info exists flags(-uniform)] } {
    set pin_aware 0
  }

  pne::set_halo_cmd $halo_x $halo_y $pin_aware
}

sta::define_cmd_args "set_pine_mp_macro_halo" \
  { -macro_name name -halo halo_list }

proc set_pine_mp_macro_halo { args } {
  sta::parse_key_args "set_pine_mp_macro_halo" args \
    keys {-macro_name -halo} \
    flags {}

  if { ![info exists keys(-macro_name)] } {
    utl::error PNE 112 "Missing required argument -macro_name"
  }

  if { ![info exists keys(-halo)] } {
    utl::error PNE 113 "Missing required argument -halo"
  }

  set macro_name $keys(-macro_name)
  set halo_list $keys(-halo)
  set n [llength $halo_list]

  if { $n == 2 } {
    # Symmetric: {horiz vert} → {horiz vert horiz vert}
    set hx [ord::microns_to_dbu [lindex $halo_list 0]]
    set vy [ord::microns_to_dbu [lindex $halo_list 1]]
    pne::set_macro_halo_cmd $macro_name $hx $vy $hx $vy
  } elseif { $n == 4 } {
    # Directional: {left bottom right top}
    set left   [ord::microns_to_dbu [lindex $halo_list 0]]
    set bottom [ord::microns_to_dbu [lindex $halo_list 1]]
    set right  [ord::microns_to_dbu [lindex $halo_list 2]]
    set top    [ord::microns_to_dbu [lindex $halo_list 3]]
    pne::set_macro_halo_cmd $macro_name $left $bottom $right $top
  } else {
    utl::error PNE 114 \
      "Halo list must have 2 values {horiz vert} or 4 values {left bottom right top}"
  }
}

namespace eval pne {
  # Internal utility functions can be added here
}

sta::define_cmd_args "set_pine_mp_soft_macros" \
  { [-utilization u] [-aspect_ratio r] }

proc set_pine_mp_soft_macros { args } {
  sta::parse_key_args "set_pine_mp_soft_macros" args \
    keys {-utilization -aspect_ratio} \
    flags {}

  # Defaults
  set utilization 0.7
  set aspect_ratio 1.0

  if { [info exists keys(-utilization)] } {
    set utilization $keys(-utilization)
  }

  if { [info exists keys(-aspect_ratio)] } {
    set aspect_ratio $keys(-aspect_ratio)
  }

  pne::set_soft_macros_cmd 1 $utilization $aspect_ratio
}

sta::define_cmd_args "pine_mp_report_soft_macros" {}

proc pine_mp_report_soft_macros { args } {
  sta::parse_key_args "pine_mp_report_soft_macros" args \
    keys {} \
    flags {}

  pne::report_soft_macros_cmd
}
