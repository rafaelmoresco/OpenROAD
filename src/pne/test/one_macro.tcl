# Basic PineMP test
# Single-cell design smoke test.
source "pine_test_utils.tcl"

read_lef "./Nangate45/Nangate45.lef"
read_def "./one_macro.def"

pne_set_quick_defaults
pne_run_and_save one_macro

puts "pass"
