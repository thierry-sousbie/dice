#ifndef __ORIENTATION_PREDICATE_BASE_2D_HXX__
#define __ORIENTATION_PREDICATE_BASE_2D_HXX__

#include "../filterType.hxx"
#include "orientationBasePrototype.hxx"

#include "../../../tools/wrappers/boostMultiprecisionFloat128.hxx"
#include "../../../tools/helpers/helpers.hxx"

#ifdef HAVE_BOOST
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#endif
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#endif

#include "../../../internal/namespace.header"

namespace internal {

  template <>
  class OrientationBaseT<2,predicate::filterType::Raw>
  {
  public:
    static const int NDIM=2;

    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r)		    
    {  
      // Ensure that we always compute in the same order ...
      int branch = (q[0]<p[0]);
      if (q[0]==p[0]) branch = (q[1]<p[1]);

      if (branch)
	{     
	  CT pqx = q[0] - p[0];
	  CT pqy = q[1] - p[1];
	  CT prx = r[0] - p[0];
	  CT pry = r[1] - p[1];

	  CT det = pqx * pry - pqy * prx;
      
	  if (det==0)
	    {
	      if (p[0]==q[0]) return (q[1]>=p[1]);
	      else return (q[0]<p[0]);

	      // if (pqx==0) return (pqy>=0);
	      // return (pqx<0);
	    }

	  return det>0;
	}
      else 
	{
	  CT pqx = p[0] - q[0];
	  CT pqy = p[1] - q[1];
	  CT prx = r[0] - q[0];
	  CT pry = r[1] - q[1];

	  CT det = pqx * pry - pqy * prx;
      
	  if (det==0)
	    {
	      if (p[0]==q[0]) return (p[1]<q[1]);
	      else return (p[0]>=q[0]);

	      // if (pqx==0) return (pqy<0);
	      // return (pqx>=0);
	    }

	  return det<0;
	}    
    }

    // NB: r is used as reference point for periodic boundary
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r,const G* geometry)		    
    { 
      T pq[2][NDIM]={{p[0],p[1]},{q[0],q[1]}};
      geometry->template checkCoordsConsistency<T,2,NDIM>(pq,r);
      return test<T,CT>(pq[0],pq[1],r);
    }
  };

#ifdef HAVE_BOOST
#ifdef HAVE_GMP

  template <>
  class OrientationBaseT<2,predicate::filterType::Exact>
  {
  public:
    static const int NDIM=2;

    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r)		    
    {           
      typedef boost::multiprecision::mpf_float mpfloat;

      mpfloat pqx = q[0]; pqx-=p[0];
      mpfloat pqy = q[1]; pqy-=p[1];
      mpfloat prx = r[0]; prx-=p[0];
      mpfloat pry = r[1]; pry-=p[1];

      mpfloat det = pqx * pry - pqy * prx;
      
      if (det==0)
	{
	  if (p[0]==q[0]) return (q[1]>=p[1]);
	  else return (q[0]<p[0]);	     
	}

      return det>0;
    }

    // NB: r is used as reference point for periodic boundary
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r,const G* geometry)		    
    {
      typedef boost::multiprecision::mpf_float mpfloat;
      mpfloat mp_pq[2][NDIM]={{p[0],p[1]},{q[0],q[1]}};
      mpfloat mp_r[2]={r[0],r[1]};
      geometry->template checkCoordsConsistency<mpfloat,2,NDIM>(mp_pq,mp_r);
      return test(mp_pq[0],mp_pq[1],mp_r);
    }
  };

#endif
#endif

  template <>
  class OrientationBaseT<2,predicate::filterType::Adaptive>
  {    
    
  public:
    static const int NDIM=2;

    template <class T, class CT=T>
    static int testAccuracy(const T p[NDIM],const T q[NDIM],
			    const T *r, CT minAcc=0)
    {  
      CT pqx = q[0] - p[0];
      CT pqy = q[1] - p[1];
      CT prx = r[0] - p[0];
      CT pry = r[1] - p[1];

      CT maxx = fabs(pqx);
      CT maxy = fabs(pqy);
      CT aprx = fabs(prx);
      CT apry = fabs(pry);

      if (maxx < aprx) maxx = aprx;
      if (maxy < apry) maxy = apry;
           
      // => See BURNIKEL & al. (95?)
      //CT eps = 1.34e-15 * maxx * maxy;
      CT eps = 3E-14 * (maxx * maxy + minAcc*minAcc); // 3.E-14
      
      CT det = pqx * pry - pqy * prx;
     
      if (det>eps) return 1;
      if (det<-eps) return 0;
      return -1;
    }

    // NB: r is used as reference point for periodic boundary
    // WARNING: this may fail for periodic boundaries if (pqr) crosses the boundary
    // AND T is long double AND r[i] is the first long double after/before p[i] or q[i]
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r,const G* geometry)		    
    { 
      T pq[2][NDIM]={{p[0],p[1]},{q[0],q[1]}};
      int res;

      // This is because we loose accuracy if coords are non consistent, so we first test
      // whether accuracy is not sufficient, we want to periodized with an mpfloat ...
      if (!geometry->template coordsAreConsistent<T,2,NDIM>(pq,r))
	{	 
	  // T should be second so that MPT_I==T when T has same precision has long double
	  // typedef typename hlp::MostPreciseFP<long double, T>
	  //   ::Result MPT_I;
	  // typedef typename hlp::IF_<hlp::SameType<T,MPT_I>::value,Float128OrMore,MPT_I>
	  //   ::Result MPT;

	  T pq2[2][NDIM]={{p[0],p[1]},{q[0],q[1]}};
	  T r2[NDIM]={r[0],r[1]};
	  
	  bool modified=geometry->template checkCoordsConsistency<T,2,NDIM>(pq2,r2);
	  
	  // If the simplex crossed a periodic boundary, we must be carefull because if 
	  // it was periodized to a coordinate very close to 0, the measured accuracy
	  // will be overestimated !
	  if (modified)
	    res=testAccuracy<T,CT>(pq2[0],pq2[1],r2,geometry->getBBoxSize());
	  else
	    res=testAccuracy<T,CT>(pq2[0],pq2[1],r2);
	  
	  //std::copy_n(&pq2[0][0],NDIM*2,&pq[0][0]);
	}
      else res=testAccuracy<T,CT>(pq[0],pq[1],r);
    
      if (res<0) 
	return OrientationBaseT<NDIM,predicate::filterType::Exact>::test(p,q,r,geometry);
      else return res;
    }
   
    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],
		    const T *r)		    
    { 
      T pq[2][NDIM]={{p[0],p[1]},{q[0],q[1]}};     
      int res=testAccuracy<T,CT>(pq[0],pq[1],r);
      
      if (res<0) 
	return OrientationBaseT<NDIM,predicate::filterType::Exact>::test(p,q,r);
      else return res;
    }
  };

}

#include "../../../internal/namespace.footer"
#endif
