#ifndef __DICE_PPROF_WRAPPER_HXX__
#define __DICE_PPROF_WRAPPER_HXX__

#ifdef USE_PPROF
#include "google/profiler.h"

#define PPROF_START_DEFAULT {ProfilerStart("/home/thierry/pprof/vps.prof");}
#define PPROF_START(x) {ProfilerStart( x );}
#define PPROF_STOP {ProfilerStop();}

#else

#define PPROF_START_DEFAULT {}
#define PPROF_START(x) {}
#define PPROF_STOP {}

#endif

#endif
