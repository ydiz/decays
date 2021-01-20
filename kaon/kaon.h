#pragma once

#include <headers/headers.h>
#include "util.h"

// #ifdef BNL
// #include "/sdcc/u/ydzhao/A2AGrid/io.h"
// #endif

#include <qlat/qlat.h>

#ifdef BNL
#include "env_bnl.h"
#endif

#ifdef CUTH_FREE_FIELD
#include "env_free_field.h"
#endif

#include "convolution.h"

