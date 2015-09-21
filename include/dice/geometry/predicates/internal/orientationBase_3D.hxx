#ifndef __ORIENTATION_PREDICATE_BASE_3D_HXX__
#define __ORIENTATION_PREDICATE_BASE_3D_HXX__

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
  class OrientationBaseT<3,predicate::filterType::Raw>
  {
  public:
    static const int NDIM=3;

    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s)		    
    {       
      CT psx = s[0] - p[0];
      CT psy = s[1] - p[1];
      CT psz = s[2] - p[2];

      CT pqx = q[0] - p[0];
      CT pqy = q[1] - p[1];
      CT pqz = q[2] - p[2];
      CT prx = r[0] - p[0];
      CT pry = r[1] - p[1];
      CT prz = r[2] - p[2];         

      CT det1 = pqy*prz - pqz*pry;
      CT det2 = pqz*prx - pqx*prz;
      CT det3 = pqx*pry - pqy*prx;
      CT det = psx*det1 + psy*det2 + psz*det3;
     
      // Simulation of simplicity ... 
      // by convention, the coordinates of p along dimension i (0<=i<NDIM) are perturbed as:
      // 
      //             Xi -> Xi - pow(eps,NDIM-i+1)
      if (det==0)
	{
	  if (det3==0) 
	    {
	      if (det2==0)
		{
		  if (det1==0)
		    {
		      // Higher order terms, check correctness ?
		      CT tmp=(psx-prx)-(psx-pqx)+(prx-pqx);
		      if (tmp==0)
			{
			  tmp=-((psy-pry)-(psy-pqy)+(pry-pqx));
			  if (tmp==0)
			    {
			      tmp=(psz-prz)-(psz-pqz)+(prz-pqz);
			      if (tmp==0) tmp=1;
			    }
			  return (tmp<0);			   
			} else return (tmp<0);
		      //
		    }
		  else return (det1<0);
		}
	      else return (det2<0);
	    }
	  else return (det3<0);
	}

      return det>0;    
    }

    // NB: s is used as reference point for periodic boundary
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s,const G* geometry)	
    {
      T pqr[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};
      geometry->template checkCoordsConsistency<T,3,NDIM>(pqr,s);
      return test<T,CT>(pqr[0],pqr[1],pqr[2],s);
    }
  }; 

#ifdef HAVE_BOOST
#ifdef HAVE_GMP

  template <>
  class OrientationBaseT<3,predicate::filterType::Exact>
  {
  public:
    static const int NDIM=3;

    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s)		    
    {  
      typedef boost::multiprecision::mpf_float mpfloat;
      
      mpfloat pqx = q[0];pqx-= p[0];
      mpfloat pqy = q[1];pqy-= p[1];
      mpfloat pqz = q[2];pqz-= p[2];

      mpfloat prx = r[0];prx-= p[0];
      mpfloat pry = r[1];pry-= p[1];
      mpfloat prz = r[2];prz-= p[2];
   
      mpfloat psx = s[0];psx-= p[0];
      mpfloat psy = s[1];psy-= p[1];
      mpfloat psz = s[2];psz-= p[2];

      mpfloat det1 = pqy*prz - pqz*pry;
      mpfloat det2 = pqz*prx - pqx*prz;
      mpfloat det3 = pqx*pry - pqy*prx;

      mpfloat det = psx*det1 + psy*det2 + psz*det3;

      // Simulation of simplicity ... 
      // by convention, the coordinates of p along dimension i (0<=i<NDIM) are perturbed as:
      // 
      //             Xi -> Xi - pow(eps,NDIM-i+1)
      if (det==0)
	{
	  if (det3==0) 
	    {
	      if (det2==0)
		{
		  if (det1==0)
		    {
		      // Higher order terms, check correctness ?
		      mpfloat tmp=(psx-prx)-(psx-pqx)+(prx-pqx);
		      if (tmp==0)
			{
			  tmp=-((psy-pry)-(psy-pqy)+(pry-pqx));
			  if (tmp==0)
			    {
			      tmp=(psz-prz)-(psz-pqz)+(prz-pqz);
			      if (tmp==0) tmp=1;
			    }
			  return (tmp<0);			   
			} else return (tmp<0);
		      //
		    }
		  else return (det1<0);
		}
	      else return (det2<0);
	    }
	  else return (det3<0);
	}

      return det>0;
    }

    // NB: s is used as reference point for periodic boundary
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s,const G* geometry)	
    {
      typedef boost::multiprecision::mpf_float mpfloat;
      mpfloat mp_pqr[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};
      mpfloat mp_s[NDIM]={s[0],s[1],s[2]};
      geometry->template checkCoordsConsistency<mpfloat,3,NDIM>(mp_pqr,mp_s);
      return test(mp_pqr[0],mp_pqr[1],mp_pqr[2],mp_s);
    }
  }; 

#endif
#endif

  template <>
  class OrientationBaseT<3,predicate::filterType::Adaptive>
  {
  public:
    static const int NDIM=3;

    template <class T, class CT=T>
    static int testAccuracy(const T p[NDIM],const T q[NDIM],const T r[NDIM],
			    const T *s, CT minAcc=0)		    
    {  
      CT psx = s[0] - p[0];
      CT psy = s[1] - p[1];
      CT psz = s[2] - p[2];

      CT pqx = q[0] - p[0];
      CT pqy = q[1] - p[1];
      CT pqz = q[2] - p[2];

      CT prx = r[0] - p[0];
      CT pry = r[1] - p[1];
      CT prz = r[2] - p[2];   

      CT maxx = fabs(pqx);
      CT maxy = fabs(pqy);
      CT maxz = fabs(pqz);

      CT aprx = fabs(prx);
      CT apsx = fabs(psx);

      CT apry = fabs(pry);
      CT apsy = fabs(psy);

      CT aprz = fabs(prz);
      CT apsz = fabs(psz);

      if (maxx < aprx) maxx = aprx;
      if (maxx < apsx) maxx = apsx;
      if (maxy < apry) maxy = apry;
      if (maxy < apsy) maxy = apsy;
      if (maxz < aprz) maxz = aprz;
      if (maxz < apsz) maxz = apsz;

      // => See BURNIKEL & al. (95?)
       // Value used in CGAL @ CGAL/internal/Static_filters/orientationBaseT_?.h 
      // is 5.1107127829973299e-15
      CT eps = 5.2e-14 * (maxx * maxy * maxz + minAcc*minAcc*minAcc);
      
      CT det1 = pqy*prz - pqz*pry;
      CT det2 = pqz*prx - pqx*prz;
      CT det3 = pqx*pry - pqy*prx;
      CT det = psx*det1 + psy*det2 + psz*det3;
      /*
      if (findFacet(0.32892722416520081552,0.26668924819391603975,1.0,
		    0.32883364621772448455,0.23839517430931708719,1.0,
		    0.296621977672874515,0.23828467614704909594,1.0,
		    p,q,r,s,2.E-5)||
	  findFacet(0.32892722416520081552,0.26668924819391603975,-1.0,
		    0.32883364621772448455,0.23839517430931708719,-1.0,
		    0.296621977672874515,0.23828467614704909594,-1.0,
		    p,q,r,s,2.E-5))
      // if (findFacet(0.32874547297384737465,0.21201735860308154602,-1.0,
      // 		    0.32883364621772448455,0.23839517430931708719,-1.0,
      // 		    0.36299584310812488264,0.2385045243585102448,-1.0,
      // 		    p,q,r,s,2.E-5)||
      // 	  findFacet(0.32874547297384737465,0.21201735860308154602,1.0,
      // 		    0.32883364621772448455,0.23839517430931708719,1.0,
      // 		    0.36299584310812488264,0.2385045243585102448,1.0,
      // 		    p,q,r,s,2.E-5))
	{
	  int val=-1;
	  if (det > eps) val=1;
	  if (det < -eps) val=0;
	  printf("Testing facet accuracy : det=%.20lg eps=%.20lg => %d\n",det,eps,val);
	  //std::cout << "Testing facet accuracy : det=" << det << "" 
	  //printf("Testing facet : det")
	}
	if ((det>-eps)&&(det<eps)) 
	glb::console->printFlush<LOG_STD>("Exact predicate required det=%e (eps=%e)!\n",det,eps);
      */	 
      


      if (det > eps) return 1;
      if (det < -eps) return 0;
      return -1;
    }
    
    // NB: s is used as reference point for periodic boundary
    // WARNING: this may fail for periodic boundaries if (pqrs) crosses the boundary
    // AND T is long double AND s[i] is the first long double 
    // after/before p[i] or q[i] or r[i]
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s,const G* geometry)	
    {
      T pqr[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};
      int res;

      // This is because we loose accuracy if coords are non consistent, so we first test
      // whether accuracy is not sufficient, we want to periodized with an mpfloat ...
      if (!geometry->template coordsAreConsistent<T,3,NDIM>(pqr,s))
	{
	   // T should be second so that MPT_I==T when T has same precision has long double
	  // typedef typename hlp::MostPreciseFP<long double, T>
	  //   ::Result MPT_I;
	  // typedef typename hlp::IF_<hlp::SameType<T,MPT_I>::value,Float128OrMore,MPT_I>
	  //   ::Result MPT;

	  // MPT pqr2[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};
	  // MPT s2[NDIM]={s[0],s[1],s[2]};
	  // geometry->template checkCoordsConsistency<MPT,3,NDIM>(pqr2,s2);
	  // res=testAccuracy(pqr2[0],pqr2[1],pqr2[2],s2);
	  
	  T pqr2[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};
	  T s2[NDIM]={s[0],s[1],s[2]};
	  
	  bool modified=geometry->template checkCoordsConsistency<T,3,NDIM>(pqr2,s2);

	  // If the simplex crossed a periodic boundary, we must be carefull because if 
	  // it was periodized to a coordinate very close to 0, the measured accuracy
	  // will be overestimated !
	  if (modified)
	    res=testAccuracy<T,CT>(pqr2[0],pqr2[1],pqr2[2],s2,geometry->getBBoxSize());
	  else
	    res=testAccuracy<T,CT>(pqr2[0],pqr2[1],pqr2[2],s2);
	}
      else res=testAccuracy<T,CT>(pqr[0],pqr[1],pqr[2],s);
     
      if (res<0)
	res=OrientationBaseT<NDIM,predicate::filterType::Exact>::test(p,q,r,s,geometry);
     
      return res;
    }

    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],
		    const T *s)	
    {
      T pqr[3][NDIM]={{p[0],p[1],p[2]},{q[0],q[1],q[2]},{r[0],r[1],r[2]}};      
      int res=testAccuracy<T,CT>(pqr[0],pqr[1],pqr[2],s);
      if (res<0) 
	return OrientationBaseT<NDIM,predicate::filterType::Exact>::test(p,q,r,s);
      else return res;
    }
  private:
    template <typename T1>
    static T1 myAbs(T1 val)
    {
      if (val<0) return -val;
      return val;
    }

    static bool findVertex(double x1,double y1,double z1, 
			   double vx, double vy, double vz, 
			   double tol=2.E-5)
    {
      if ((myAbs(vx-x1)<=tol)&&
	  (myAbs(vy-y1)<=tol)&&
	  (myAbs(vz-z1)<=tol))
	{
	  return true;
	}
      return false;
    }

    template <typename T1, typename T2>
    static bool findFacet(T1 x1,T1 y1,T1 z1,
			  T1 x2,T1 y2,T1 z2, 
			  T1 x3,T1 y3,T1 z3,
			  T2 p[NDIM], T2 q[NDIM], T2 r[NDIM], T2 s[NDIM], 
			  double tol=2.E-5)
    {
      int nFound=0;
      bool found=false;
      found|=findVertex(x1,y1,z1,p[0],p[1],p[2],tol);
      found|=findVertex(x1,y1,z1,q[0],q[1],q[2],tol);
      found|=findVertex(x1,y1,z1,r[0],r[1],r[2],tol);
      found|=findVertex(x1,y1,z1,s[0],s[1],s[2],tol);
      if (found) nFound++;
      found=false;
      found|=findVertex(x2,y2,z2,p[0],p[1],p[2],tol);
      found|=findVertex(x2,y2,z2,q[0],q[1],q[2],tol);
      found|=findVertex(x2,y2,z2,r[0],r[1],r[2],tol);
      found|=findVertex(x2,y2,z2,s[0],s[1],s[2],tol);
      if (found) nFound++;
      found=false;
      found|=findVertex(x3,y3,z3,p[0],p[1],p[2],tol);
      found|=findVertex(x3,y3,z3,q[0],q[1],q[2],tol);
      found|=findVertex(x3,y3,z3,r[0],r[1],r[2],tol);
      found|=findVertex(x3,y3,z3,s[0],s[1],s[2],tol);
      if (found) nFound++;

      return (nFound==3);      
    }

  }; 

}

#include "../../../internal/namespace.footer"
#endif
