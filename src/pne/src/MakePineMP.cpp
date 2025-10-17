// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/MakePineMP.h"

#include <tcl.h>

#include "utl/decode.h"

extern "C" {
extern int Pne_Init(Tcl_Interp* interp);
}

namespace pne {

extern const char* pne_tcl_inits[];

void initPineMP(Tcl_Interp* tcl_interp)
{
  Pne_Init(tcl_interp);
  utl::evalTclInit(tcl_interp, pne::pne_tcl_inits);
}

}  // namespace pne
