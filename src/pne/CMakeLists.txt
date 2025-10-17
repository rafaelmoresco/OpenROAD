# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

include("openroad")

add_library(pne_lib
  src/PineMP.cpp
)

target_include_directories(pne_lib
  PUBLIC
    include
)

target_link_libraries(pne_lib
    utl_lib
    dbSta_lib
    OpenSTA
    odb
)

swig_lib(NAME      pne
         NAMESPACE pne
         I_FILE    src/pne.i
         SCRIPTS   src/pne.tcl
)

target_sources(pne
  PRIVATE
    src/MakePineMP.cpp
)

messages(
  TARGET pne
)

target_include_directories(pne
  PUBLIC
    include
  PRIVATE
    src
)

target_link_libraries(pne
    pne_lib
    gui
)

if(ENABLE_TESTS)
  add_subdirectory(test)
endif()