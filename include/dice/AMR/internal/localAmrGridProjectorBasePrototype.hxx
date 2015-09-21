#ifndef __LOCAL_AMR_GRID_PROJECTOR_BASE_PROTOTYPE_HXX__
#define __LOCAL_AMR_GRID_PROJECTOR_BASE_PROTOTYPE_HXX__

#include "../../dice_globals.hxx"

#include "../../internal/namespace.header"

namespace internal {
  template <int ND, class AMR, class MESH, 
	    typename F, typename HF, 
	    bool checkAccuracy, 
	    int dummy=0>
  class LocalAmrGridProjectorBaseT;
}

#include "../../internal/namespace.footer"
#endif
