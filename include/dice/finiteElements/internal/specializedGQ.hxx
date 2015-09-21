#ifndef __SPECIALIZED_GAUSS_QUADRATURE_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_HXX__

#include "../finiteElementTypeEnum.hxx"

#include "../../internal/namespace.header"

namespace internal {
  namespace gaussQuadrature {
   
    template <int ND, int DEG, fETypeE::Type FET>
    class ERROR_QUADRATURE_RULE_NOT_YET_IMPLEMENTED;

    // Prototype
    template <int ND, int DEG, fETypeE::Type FET>
    struct SpecializedGQT
    {
      static const int NP  = 0;

      template <class F>
      static double get(F& functor)
      {
	ERROR_QUADRATURE_RULE_NOT_YET_IMPLEMENTED<ND,DEG,FET> error;
	return 0;
      }
    };
  
  }
}

#include "../../internal/namespace.footer"

#include "specializedGQ_simplex.hxx"
#include "specializedGQ_cuboid.hxx"

#endif
