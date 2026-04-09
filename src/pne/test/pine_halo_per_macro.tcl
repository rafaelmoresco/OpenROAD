# Test per-macro directional halo override.
source "pine_test_utils.tcl"

pne_load_design "macro_only.def" "macro_only.lef"

# Global halo baseline.
set_pine_mp_halo -halo_x 3 -halo_y 3

# Per-macro override with asymmetric directional halo:
# {left bottom right top} in microns.
set_pine_mp_macro_halo -macro_name U1 -halo {8 4 2 1}

# Symmetric per-macro override: {horiz vert}.
set_pine_mp_macro_halo -macro_name U6 -halo {6 3}

pne_set_quick_defaults
pne_run_and_save pine_halo_per_macro

pne_assert_macros_inside_core

puts "pass"
