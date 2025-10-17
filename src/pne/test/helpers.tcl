# Helper functions for PineMP tests

if { [info exists ::env(TEST_TMPDIR)] } {
  set test_dir $::env(TEST_TMPDIR)
} else {
  set test_dir [file dirname [file normalize [info script]]]
}

if { [info exists ::env(RESULTS_DIR)] } {
  set result_dir $::env(RESULTS_DIR)
} else {
  set result_dir [file join $test_dir "results"]
}

proc make_result_dir { } {
  variable result_dir
  if { ![file exists $result_dir] } {
    file mkdir $result_dir
  }
  return $result_dir
}

proc make_result_file { filename } {
  variable result_dir

  make_result_dir

  set root [file rootname $filename]
  set ext [file extension $filename]
  set filename "$root-tcl$ext"
  return [file join $result_dir $filename]
}

proc diff_file { expected_file actual_file } {
  if { [file exists $expected_file] && [file exists $actual_file] } {
    if { [catch { exec diff $expected_file $actual_file } result] } {
      puts "Files differ:"
      puts $result
      return 1
    } else {
      puts "Files match"
      return 0
    }
  } else {
    puts "One or both files missing"
    return 1
  }
}

# suppress common messages
suppress_message ODB 127
suppress_message ODB 134
suppress_message ORD 30
