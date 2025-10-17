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

namespace eval pne {
  # Internal utility functions can be added here
}
