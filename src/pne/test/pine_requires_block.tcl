source "pine_test_utils.tcl"

pne_expect_error \
  "pine_mp requires a loaded block" \
  { pine_mp } \
  "*No block found for PineMP placement*"

puts "pass"
