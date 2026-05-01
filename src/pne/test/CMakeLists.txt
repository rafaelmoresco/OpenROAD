# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021-2025, The OpenROAD Authors

or_integration_tests(
  "pne"
  PASSFAIL_TESTS
    one_macro
    pine_basic
    macro_only
    no_unfixed_macros
    guides1
    fixed_ios1
    io_constraints1
    io_constraints6
    mixed_ios1
    io_pads1
    orientation_improve1
    boundary_push1
    boundary_push2
    centralization1
    clocked_macro
    fixed_macros
    pine_api_defaults
    pine_api_required_args
    pine_api_strategies
    pine_repeat_run
    pine_requires_block
    pine_halo_basic
    pine_halo_uniform
    pine_halo_per_macro
    pine_soft_macros_api
)
