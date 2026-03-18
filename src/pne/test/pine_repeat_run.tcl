source "pine_test_utils.tcl"

pne_load_design "boundary_push1.def" "orientation_improve1.lef"
pne_set_quick_defaults

if { ![pine_mp -num_threads 1] } {
  error "first pine_mp invocation failed"
}
if { ![pine_mp -num_threads 1] } {
  error "second pine_mp invocation failed"
}

set def_file [make_result_file "pine_repeat_run.def"]
write_def $def_file
if { ![file exists $def_file] } {
  error "missing output DEF for repeated run"
}

puts "pass"
