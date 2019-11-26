#pragma once


#define PARALLEL_FOR_LOOP        _Pragma("omp parallel for schedule(static)")
#define parallel_for       PARALLEL_FOR_LOOP for

#include "io.h"
#include "constants_macro.h"
#include "pGG.h"
#include "utils.h"
#include "../qlat_wrapper/qlat_wrapper.h"


