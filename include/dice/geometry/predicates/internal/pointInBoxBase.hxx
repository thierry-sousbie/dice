#ifndef __POINT_IN_BOX_BASE_HXX__
#define __POINT_IN_BOX_BASE_HXX__

#include "../../geometricProperties.hxx"

#include "../../../tools/wrappers/boostMultiprecisionFloat128.hxx"

#ifdef HAVE_BOOST
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#endif
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#endif


#include "../../../internal/namespace.header"

namespace internal {
  template <int ND, int filterType>
  class PointInBoxBaseT
  {
  public:
    
    template <class T>
    static int test(const T box[][ND], const T pCoord[ND])
    {
      for (int i=0;i<ND;++i)
	{
	  if (box[0][i]<box[1][i])
	    {
	      if ((pCoord[i]<=box[0][i])||(pCoord[i]>box[1][i])) return 0;
	    }
	  else
	    {
	      if ((pCoord[i]<=box[1][i])||(pCoord[i]>box[0][i])) return 0;
	    }
	}
      return 1;
    }

    template <class T, class G>
    static int test(const T box[][ND], const T pCoord[ND], const G* geometry)
    {
      if (!geometry->template coordsAreConsistent<T,2,ND>(box,pCoord))
	{
	  T pq2[2][ND];
	  T r2[ND];

	  for (int i=0;i<ND;++i)
	    r2[i]=pCoord[i];
	  for (int i=0;i<ND;++i) 
	    {
	      pq2[0][i]=box[0][i];
	      pq2[1][i]=box[1][i];
	      r2[i]=pCoord[i];
	    }
	  
	  geometry->template checkCoordsConsistency<T,2,ND>(pq2,r2);
	    
	  if (filterType==predicate::filterType::Raw)
	    return test(pq2,r2);

	  for (int i=0;i<ND;++i)
	    {
	      T eps=1.E-15*geometry->getBBoxSize(i);
	      if ((fabs(pq2[0][i]-r2[i])<eps)||
		  (fabs(pq2[1][i]-r2[i])<eps))
		{
		  // We need exact computations !
		  typedef boost::multiprecision::mpf_float mpfloat;
		  mpfloat mp_pq[2][ND];
		  mpfloat mp_r[ND];
		  
		  for (int i=0;i<ND;++i)
		    mp_r[i]=pCoord[i];
		  for (int i=0;i<ND;++i) 
		    {
		      mp_pq[0][i]=box[0][i];
		      mp_pq[1][i]=box[1][i];
		      mp_r[i]=pCoord[i];
		    }
		  geometry->template checkCoordsConsistency<mpfloat,2,ND>(mp_pq,mp_r);
		  return test(mp_pq,mp_r);
		}
	    }

	  return test(pq2,r2);
	}
      else return test(box,pCoord);
    }
    
  };

} // namespace internal
 
#include "../../../internal/namespace.footer"
#endif
