# Basic PineMP test
source "helpers.tcl"
read_lef Nangate45/Nangate45.lef
read_def basic.def

# Test basic PineMP command
pine_mp -num_threads 1

# Test debug command
pine_mp_debug -debug

# Run placement again with debug enabled
pine_mp -num_threads 1

set def_file [make_result_file basic.def]
write_def $def_file
diff_file basic.defok $def_file
