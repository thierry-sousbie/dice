#ifndef __STATIC_POTENTIAL_COLDICE_CONFIG_H__
#define __STATIC_POTENTIAL_COLDICE_CONFIG_H__

#include "../config.h"

// boundary conditions => boxed
#ifdef D_BOUNDARY_TYPE
#undef D_BOUNDARY_TYPE
#define D_BOUNDARY_TYPE BOXED
#endif

#endif
