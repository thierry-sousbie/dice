#ifndef __BOOST_MULTIPRECISION_FLOAT128_FIX_HXX__
#define __BOOST_MULTIPRECISION_FLOAT128_FIX_HXX__

/**
 * @file 
 * @brief This header should be included in place of <boost/multiprecision/float128.hpp> in
 * order to fix a bug in boost that prevents casting from float128 to mp_float 
 * type (tested for Boost v1.57.0). It also a defines a Float128OrMore type that is 
 * available also when libquadmath is not used and is guaranted to have at least 
 * 128bit precision.
 * @author Thierry Sousbie
 */

#include <boost/cstdfloat.hpp>
#include "../helpers/helpers.hxx"
#ifdef HAVE_QUADMATH

#include <boost/multiprecision/float128.hpp>

namespace boost {
  namespace multiprecision {
    namespace backends {
      
      inline bool eval_is_zero(const float128_backend& val) BOOST_NOEXCEPT
      {      
	return boost::multiprecision::default_ops::eval_is_zero(val);
      }

      inline int eval_get_sign(const float128_backend& val) BOOST_NOEXCEPT
      {      
	return boost::multiprecision::default_ops::eval_get_sign(val);
      }
           
    }
  }
}

#include "../../internal/namespace.header"
typedef boost::multiprecision::float128 Float128OrMore;
typedef hlp::IsTrue Float128OrMore_IsQuadruplePrecision;
#include "../../internal/namespace.footer"

#else // NO quadmath

#include <boost/multiprecision/gmp.hpp>

#include "../../internal/namespace.header"
typedef boost::multiprecision::number<boost::multiprecision::gmp_float<35> > Float128OrMore;//mpf_float_50
typedef hlp::IsFalse Float128OrMore_IsQuadruplePrecision;
#include "../../internal/namespace.footer"

#endif

#endif
