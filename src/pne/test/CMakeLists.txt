# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021-2025, The OpenROAD Authors

add_executable(pne_test
  regression
)

target_link_libraries(pne_test
  pne
  utl_lib
)