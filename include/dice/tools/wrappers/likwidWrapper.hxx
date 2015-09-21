#ifndef __DICE_LIKWID_WRAPPER_HXX__
#define __DICE_LIKWID_WRAPPER_HXX__

#ifdef USE_LIKWID
#include "likwid.h"
#else

#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(reg)
#define LIKWID_MARKER_STOP(reg)
#define LIKWID_MARKER_CLOSE

#endif

#endif
