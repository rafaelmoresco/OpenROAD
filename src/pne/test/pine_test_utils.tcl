# Helper functions for PineMP pass/fail regressions.

source "helpers.tcl"

proc pne_expect_error { description script expected_pattern } {
  set status [catch { uplevel 1 $script } message]
  if { $status == 0 } {
    error "$description: expected an error"
  }
  if { ![string match $expected_pattern $message] } {
    error "$description: expected '$expected_pattern', got '$message'"
  }
}

proc pne_load_design { def_name { macro_lef "" } { io_lef "" } } {
  read_lef "./Nangate45/Nangate45.lef"

  if { $macro_lef != "" } {
    read_lef [file join "." "testcases" $macro_lef]
  }

  if { $io_lef != "" } {
    read_lef [file join "." "Nangate45_io" $io_lef]
  }

  read_def [file join "." "testcases" $def_name]
}

proc pne_set_quick_defaults { { pin_strategy "connectivity" } } {
  # Keep runtime bounded while still exercising iterative optimization.
  set_pine_mp_iterations 2
  set_pine_mp_initial_weights -internal_weight 0.8 -io_weight 0.2
  set_pine_mp_final_weights -internal_weight 0.5 -io_weight 0.5
  set_pine_mp_sa_params -initial_temp 250.0 -cooling_rate 0.9 -max_iterations 200
  set_pine_mp_pin_strategy $pin_strategy
}

proc pne_run_and_save { test_name { num_threads 1 } } {
  if { ![pine_mp -num_threads $num_threads] } {
    error "pine_mp returned false for $test_name"
  }

  set def_file [make_result_file "$test_name.def"]
  write_def $def_file
  if { ![file exists $def_file] } {
    error "missing output DEF for $test_name"
  }

  return $def_file
}

proc pne_run_scenario { test_name def_name { macro_lef "" } { io_lef "" } } {
  pne_load_design $def_name $macro_lef $io_lef
  pne_set_quick_defaults
  pne_run_and_save $test_name
}

proc pne_assert_macros_clear_io_pads {} {
  set block [ord::get_db_block]

  set io_obstacles {}
  foreach inst [$block getInsts] {
    set master [$inst getMaster]
    if { !([$inst isPad] || [$master isCover]) } {
      continue
    }
    if { !([$inst isFixed] || [$inst isPlaced]) } {
      continue
    }

    set bbox [$inst getBBox]
    lappend io_obstacles [list [$inst getName] [$bbox xMin] [$bbox yMin] [$bbox xMax] [$bbox yMax]]
  }

  foreach inst [$block getInsts] {
    set master [$inst getMaster]
    if { !([$inst isBlock] && ![$inst isFixed] && ![$master isPad] && ![$master isCover]) } {
      continue
    }

    set bbox [$inst getBBox]
    set x0 [$bbox xMin]
    set y0 [$bbox yMin]
    set x1 [$bbox xMax]
    set y1 [$bbox yMax]

    foreach obstacle $io_obstacles {
      lassign $obstacle obstacle_name ox0 oy0 ox1 oy1
      set overlap_x [expr {max(0, min($x1, $ox1) - max($x0, $ox0))}]
      set overlap_y [expr {max(0, min($y1, $oy1) - max($y0, $oy0))}]
      if { $overlap_x > 0 && $overlap_y > 0 } {
        error "Macro [$inst getName] overlaps IO obstacle $obstacle_name"
      }
    }
  }
}

proc pne_assert_macros_inside_core {} {
  set block [ord::get_db_block]
  set core [$block getCoreArea]
  set cx0 [$core xMin]
  set cy0 [$core yMin]
  set cx1 [$core xMax]
  set cy1 [$core yMax]

  foreach inst [$block getInsts] {
    set master [$inst getMaster]
    if { !([$inst isBlock] && ![$inst isFixed] && ![$master isPad] && ![$master isCover]) } {
      continue
    }

    set bbox [$inst getBBox]
    set x0 [$bbox xMin]
    set y0 [$bbox yMin]
    set x1 [$bbox xMax]
    set y1 [$bbox yMax]

    if { $x0 < $cx0 || $y0 < $cy0 || $x1 > $cx1 || $y1 > $cy1 } {
      error "Macro [$inst getName] bbox ($x0,$y0)-($x1,$y1) exceeds core area ($cx0,$cy0)-($cx1,$cy1)"
    }
  }
}

proc pne_assert_def_origins_within_die { def_file } {
  set stream [open $def_file r]

  set die_max_x -1
  set die_max_y -1
  set in_components 0

  while { [gets $stream line] >= 0 } {
    if { [regexp {^DIEAREA\s+\(\s*[-0-9]+\s+[-0-9]+\s*\)\s+\(\s*([-0-9]+)\s+([-0-9]+)\s*\)\s*;} $line -> max_x max_y] } {
      set die_max_x $max_x
      set die_max_y $max_y
      continue
    }

    if { [regexp {^COMPONENTS\s+} $line] } {
      set in_components 1
      continue
    }
    if { [regexp {^END COMPONENTS} $line] } {
      set in_components 0
      continue
    }

    if { $in_components && [regexp {^\s*-\s+(\S+)\s+\S+\s+\+\s+PLACED\s+\(\s*([-0-9]+)\s+([-0-9]+)\s*\)} $line -> inst_name x y] } {
      if { $x < 0 || $y < 0 || $x > $die_max_x || $y > $die_max_y } {
        close $stream
        error "Instance $inst_name origin ($x,$y) exceeds die bounds (0,0)-($die_max_x,$die_max_y)"
      }
    }
  }

  close $stream

  if { $die_max_x < 0 || $die_max_y < 0 } {
    error "Could not parse DIEAREA from $def_file"
  }
}
