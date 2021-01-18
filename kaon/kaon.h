#pragma once

#include <headers/headers.h>
#include "util.h"

// #ifdef BNL
// #include "/sdcc/u/ydzhao/A2AGrid/io.h"
// #endif

#include <qlat/qlat.h>

#ifdef CUTH
#include "env_cuth.h"
#endif


#ifdef ARGONNE
#include "env_argonne.h"
#endif

#ifdef BNL
#include "env_bnl.h"
#endif

// #include "kaon_init.h"

#include "convolution.h"

