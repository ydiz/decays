#pragma once

#include <headers/headers.h>
#include "util.h"

#ifdef CORI
#include "env_cori.h"
#endif

#ifdef BNL
#include <qlat/qlat.h>
#include "env_bnl.h"
#endif

#ifdef CUTH_FREE_FIELD
#include "env_free_field.h"
#endif

#include "convolution.h"

