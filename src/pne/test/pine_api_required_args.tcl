source "pine_test_utils.tcl"

pne_expect_error \
  "initial weights require -internal_weight" \
  { set_pine_mp_initial_weights -io_weight 0.2 } \
  "*Missing required argument -internal_weight*"

pne_expect_error \
  "initial weights require -io_weight" \
  { set_pine_mp_initial_weights -internal_weight 0.8 } \
  "*Missing required argument -io_weight*"

pne_expect_error \
  "final weights require -internal_weight" \
  { set_pine_mp_final_weights -io_weight 0.5 } \
  "*Missing required argument -internal_weight*"

pne_expect_error \
  "final weights require -io_weight" \
  { set_pine_mp_final_weights -internal_weight 0.5 } \
  "*Missing required argument -io_weight*"

pne_expect_error \
  "objective requires -style" \
  { set_pine_mp_objective -spread_weight 1.0 } \
  "*Missing required argument -style*"

pne_expect_error \
  "objective style must be min_wl or spread" \
  { set_pine_mp_objective -style max_wl } \
  "*-style must be min_wl or spread*"

pne_expect_error \
  "spread weight must be non-negative" \
  { set_pine_mp_objective -style spread -spread_weight -1 } \
  "*-spread_weight*is not a positive float*"

pne_expect_error \
  "min channel takes one or two values" \
  { set_pine_mp_halo -halo_x 5 -halo_y 5 -min_channel {1 2 3} } \
  "*-min_channel must be*"

pne_expect_error \
  "min channel must be non-negative" \
  { set_pine_mp_halo -halo_x 5 -halo_y 5 -min_channel {-1 2} } \
  "*-min_channel values must be non-negative*"

puts "pass"
