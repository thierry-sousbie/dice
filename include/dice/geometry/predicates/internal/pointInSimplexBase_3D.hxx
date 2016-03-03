#ifndef __POINT_IN_SIMPLEX_BASE_3D_HXX__
#define __POINT_IN_SIMPLEX_BASE_3D_HXX__

#include <string.h> // for memcpy

#include "pointInSimplexBasePrototype.hxx"
#include "pointInSimplexBase_2D.hxx"

#include "../orientation.hxx"

#include "../../simplexInterpolator.hxx"
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

//#define DEBUG_ME 1

namespace internal {  
  
  template <int filterType>
  class PointInSimplexBaseT<3,filterType>
  {

  public:
    static const int NDIM=3;
    typedef predicate::OrientationT<3,filterType> Orientation;
    /*
    template <class T, class T2, class G, class CT=T2>
    static int getRayFacetIntersection(T vCoord[][NDIM], const T pCoord[NDIM], 
				       const int dim, T2 &zCoord, const G*geometry)
    {     
      T vCoord2D[NDIM][2];
      T pCoord2D[2];    
      
      int j=0;
      for (int i=0;i<NDIM;++i) 
	{
	  if (i!=dim) 
	    {	   
	      pCoord2D[j]=pCoord[i];
	      vCoord2D[0][j]=vCoord[0][i];
	      vCoord2D[1][j]=vCoord[1][i];
	      vCoord2D[2][j]=vCoord[2][i];
	      ++j;
	    }	 
	}

      if (test2D<T>(vCoord2D,pCoord2D,dim,geometry))
	{

	  // Sort the two faces so that the points are always in the same order, whatever 
	  // the original order of the simplices. This ensure IEEE FP errors consistency !
	  // FIXME: check the necessity of it, this is SLOW !
	  
	  for (int j=0;j<NDIM-1;++j)
	    for (int k=j+1;k<NDIM;++k)
	      {
		if (vCoord[j][0]<vCoord[k][0])
		  for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		else if (vCoord[j][0]==vCoord[k][0])
		  {
		    if (vCoord[j][1]<vCoord[k][1])
		      for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);

		    else if (vCoord[j][1]==vCoord[k][1])
		      {
			if (vCoord[j][2]<vCoord[k][2])
			  for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		      }
		  }
	      }
	  zCoord = interpolate<T,G,T2,CT>(vCoord,pCoord,dim,geometry);
	  return 1;
	}      
      return 0;
    }

    template <class T, class T2, class CT=T2>
    static int getRayFacetIntersection(T vCoord[][NDIM], const T pCoord[NDIM], 
				       const int dim, T2 &zCoord)
    {     
      T vCoord2D[NDIM][2];
      T pCoord2D[2];    
      
      int j=0;
      for (int i=0;i<NDIM;++i) 
	{
	  if (i!=dim) 
	    {	   
	      pCoord2D[j]=pCoord[i];
	      vCoord2D[0][j]=vCoord[0][i];
	      vCoord2D[1][j]=vCoord[1][i];
	      vCoord2D[2][j]=vCoord[2][i];
	      ++j;
	    }	 
	}

      if (test2D<T>(vCoord2D,pCoord2D))
	{

	  // Sort the two faces so that the points are always in the same order, whatever 
	  // the original order of the simplices. This ensure IEEE FP errors consistency !
	  // FIXME: check the necessity of it, this is SLOW !
	  
	  for (int j=0;j<NDIM-1;++j)
	    for (int k=j+1;k<NDIM;++k)
	      {
		if (vCoord[j][0]<vCoord[k][0])
		  for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		else if (vCoord[j][0]==vCoord[k][0])
		  {
		    if (vCoord[j][1]<vCoord[k][1])
		      for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);

		    else if (vCoord[j][1]==vCoord[k][1])
		      {
			if (vCoord[j][2]<vCoord[k][2])
			  for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		      }
		  }
	      }
	  
	  zCoord = interpolate<T,T2,CT>(vCoord,pCoord,dim);
	  return 1;
	}      
      return 0;
    }
    */
    // FIXME: IT IS UGLY TO DUPLICATE THIS FUNCTION, FIND A BETTER WAY ?
    template <class T, class T2, class CT=T2, bool INTERP=true>
    static int getSegmentFacetIntersection(const T vCoord[NDIM][NDIM], 
					   const T pCoord[NDIM], 
					   int dim, T otherCoord, 
					   T2 &intersectionCoord)
    {
      T vCoord2D[NDIM][2];
      T pCoord2D[2];    
      
      if (dim==0)
	{
	  pCoord2D[0]=pCoord[1];pCoord2D[1]=pCoord[2];
	  vCoord2D[0][0]=vCoord[0][1];vCoord2D[0][1]=vCoord[0][2];
	  vCoord2D[1][0]=vCoord[1][1];vCoord2D[1][1]=vCoord[1][2];
	  vCoord2D[2][0]=vCoord[2][1];vCoord2D[2][1]=vCoord[2][2];
	}
      else if (dim==1)
	{
	  pCoord2D[0]=pCoord[0];pCoord2D[1]=pCoord[2];
	  vCoord2D[0][0]=vCoord[0][0];vCoord2D[0][1]=vCoord[0][2];
	  vCoord2D[1][0]=vCoord[1][0];vCoord2D[1][1]=vCoord[1][2];
	  vCoord2D[2][0]=vCoord[2][0];vCoord2D[2][1]=vCoord[2][2];
	}
      else
	{
	  pCoord2D[0]=pCoord[0];pCoord2D[1]=pCoord[1];
	  vCoord2D[0][0]=vCoord[0][0];vCoord2D[0][1]=vCoord[0][1];
	  vCoord2D[1][0]=vCoord[1][0];vCoord2D[1][1]=vCoord[1][1];
	  vCoord2D[2][0]=vCoord[2][0];vCoord2D[2][1]=vCoord[2][1];
	}
         
      if (test2D<T>(vCoord2D,pCoord2D))
	{
	  double testPointCoord[NDIM];
	  std::copy(pCoord,pCoord+NDIM,testPointCoord);
	  int res1 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
				       testPointCoord);
      
	  testPointCoord[dim]=otherCoord;
	  int res2 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
				       testPointCoord);
	    
	  if (res1!=res2)
	    {
	      // Sort the vertices so that they are always in the same order. 
	      // This is to ensure that the interpolated intersection coordinates 
	      // are not affected by rounding errors when changing the facet's 
	      // vertices order.
	      /*
	      for (int j=0;j<NDIM-1;++j)
		for (int k=j+1;k<NDIM;++k)
		  {
		    if (vCoord[j][0]<vCoord[k][0])
		      for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		    else if (vCoord[j][0]==vCoord[k][0])
		      {
			if (vCoord[j][1]<vCoord[k][1])
			  for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);

			else if (vCoord[j][1]==vCoord[k][1])
			  {
			    if (vCoord[j][2]<vCoord[k][2])
			      for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
			  }
		      }
		  }
	      */
	      if (INTERP) intersectionCoord = interpolate<T,T2,CT>(vCoord,pCoord,dim);		
	      return 1;
	    }
	}
	
      return 0;
    }

    template <class T, class T2, class G, class CT=T2, bool INTERP=true>
    static int getSegmentFacetIntersection(const T vCoord[NDIM][NDIM], 
					   const T pCoord[NDIM], 
					   int dim, T otherCoord, 
					   T2 &intersectionCoord,
					   const G* geometry)
    {
      T vCoord2D[NDIM][2];
      T pCoord2D[2];    
      
      if (dim==0)
	{
	  pCoord2D[0]=pCoord[1];pCoord2D[1]=pCoord[2];
	  vCoord2D[0][0]=vCoord[0][1];vCoord2D[0][1]=vCoord[0][2];
	  vCoord2D[1][0]=vCoord[1][1];vCoord2D[1][1]=vCoord[1][2];
	  vCoord2D[2][0]=vCoord[2][1];vCoord2D[2][1]=vCoord[2][2];
	}
      else if (dim==1)
	{
	  pCoord2D[0]=pCoord[0];pCoord2D[1]=pCoord[2];
	  vCoord2D[0][0]=vCoord[0][0];vCoord2D[0][1]=vCoord[0][2];
	  vCoord2D[1][0]=vCoord[1][0];vCoord2D[1][1]=vCoord[1][2];
	  vCoord2D[2][0]=vCoord[2][0];vCoord2D[2][1]=vCoord[2][2];
	}
      else
	{
	  pCoord2D[0]=pCoord[0];pCoord2D[1]=pCoord[1];
	  vCoord2D[0][0]=vCoord[0][0];vCoord2D[0][1]=vCoord[0][1];
	  vCoord2D[1][0]=vCoord[1][0];vCoord2D[1][1]=vCoord[1][1];
	  vCoord2D[2][0]=vCoord[2][0];vCoord2D[2][1]=vCoord[2][1];
	}
      
      if (test2D<T>(vCoord2D,pCoord2D,dim,geometry))
	  {
	    double testPointCoord[NDIM];
	    std::copy(pCoord,pCoord+NDIM,testPointCoord);
	    int res1 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
					 testPointCoord,geometry);
      
	    testPointCoord[dim]=otherCoord;
	    int res2 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
					 testPointCoord,geometry);

	    if (res1!=res2)
	      {
		// Sort the vertices so that they are always in the same order. 
		// This is to ensure that the interpolated intersection coordinates 
		// are not affected by rounding errors when changing the facet's 
		// vertices order.
		/*
		for (int j=0;j<NDIM-1;++j)
		  for (int k=j+1;k<NDIM;++k)
		    {
		      if (vCoord[j][0]<vCoord[k][0])
			for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
		      else if (vCoord[j][0]==vCoord[k][0])
			{
			  if (vCoord[j][1]<vCoord[k][1])
			    for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);

			  else if (vCoord[j][1]==vCoord[k][1])
			    {
			      if (vCoord[j][2]<vCoord[k][2])
				for (int l=0;l<NDIM;++l) std::swap(vCoord[j][l],vCoord[k][l]);
			    }
			}
		    }
		*/
		if (INTERP) intersectionCoord = interpolate<T,G,T2,CT>(vCoord,pCoord,dim,geometry);
		return 1;
	      }
	  }
	
      return 0;
    }

    template <class T>
    static int test(const T vCoord[][NDIM], const T pCoord[NDIM])
    {
      int res=0;
      res+=Orientation::test(vCoord[0],vCoord[1],vCoord[2],pCoord);
      res+=Orientation::test(vCoord[0],vCoord[2],vCoord[3],pCoord);
      res+=Orientation::test(vCoord[0],vCoord[3],vCoord[1],pCoord);
      res+=Orientation::test(vCoord[1],vCoord[3],vCoord[2],pCoord);
      return ((res==4)||(res==0));
    }

    template <class T, class G>
    static int test(const T vCoord[][NDIM], const T pCoord[NDIM], const G* geometry)
    {
      int res=0;
      res+=Orientation::test(vCoord[0],vCoord[1],vCoord[2],pCoord,geometry);
      res+=Orientation::test(vCoord[0],vCoord[2],vCoord[3],pCoord,geometry);
      res+=Orientation::test(vCoord[0],vCoord[3],vCoord[1],pCoord,geometry);
      res+=Orientation::test(vCoord[1],vCoord[3],vCoord[2],pCoord,geometry);
      return ((res==4)||(res==0));      
    }

  private:       

    template <class T>
    static int test2D(const T vCoord[][2], const T pCoord[2])
    {
      return PointInSimplexBaseT<2,filterType>::template test<T>(vCoord,pCoord);
    }

    template <class T, class G>
    static int test2D(const T vCoord[][2], const T pCoord[2], int dim, const G* geometry)
    {
      //return PointInSimplexBaseT<2,filterType>::template test<T>(vCoord,pCoord,geometry);
      
      typedef GeometricPropertiesT<typename G::Coord,G::NDIM-1,G::NDIM_W-1,
				   G::BOUNDARY_TYPE,G::WORLD_BOUNDARY_TYPE> G2D;	
      G2D geometry2D(*geometry,dim);
      return PointInSimplexBaseT<2,filterType>::template test<T>(vCoord,pCoord,&geometry2D);      
    }
    
    template <class T, class RT=T, class CT=T>
    static RT interpolateZ(const T vCoord[][NDIM], const T pCoord[NDIM])
    {
      //static const double tol=1.e5*std::numeric_limits<CT>::epsilon();
      static const double tol=1.e4*std::numeric_limits<double>::epsilon();
      CT det=(vCoord[0][0]-vCoord[2][0]);
      det*=(vCoord[1][1]-vCoord[2][1]);
      CT det2=(vCoord[0][1]-vCoord[2][1]);
      det2*=(vCoord[1][0]-vCoord[2][0]);  

      // Check if the determinant is precise enough 
#if defined(HAVE_BOOST) && defined(HAVE_GMP)
      CT det_inv=det-det2;
      if ((fabs(det_inv)<std::max(fabs(det2),fabs(det))*tol)||(det_inv==0))
	{
	  typedef boost::multiprecision::mpf_float mpfloat;
	   mpfloat a=vCoord[0][0];a-=vCoord[2][0];
	   mpfloat b=vCoord[1][1];b-=vCoord[2][1];	   
	   a*=b;

	   b=vCoord[0][1];b-=vCoord[2][1];	   
	   mpfloat c=vCoord[1][0];c-=vCoord[2][0];	   
	   b*=c;
	   
	   // because of SOS, a degenerate face cannot be intersected, so a!=b
	   c=1;
	   if (a!=b) c/=(a-b); 
	   det_inv = c.convert_to<CT>();	  
	}
      else
	{
	  det=det_inv;
	  det2=1;
	  det_inv=det2/det;	  	  	  
	}
#else
      if (det==det2) 
	return hlp::numericStaticCast<RT>((vCoord[0][2]+vCoord[1][2]+vCoord[2][2])/3);
      det=det_inv;
      det2=1;
      det_inv=det2/det;
#endif
      
      det=(pCoord[0]-vCoord[2][0]);
      det*=(vCoord[1][1]-vCoord[2][1]);
      det2=(pCoord[1]-vCoord[2][1]);
      det2*=(vCoord[2][0]-vCoord[1][0]);
      CT L1=det_inv*(det+det2);

      det=(pCoord[0]-vCoord[2][0]);
      det*=(vCoord[2][1]-vCoord[0][1]);
      det2=(pCoord[1]-vCoord[2][1]);
      det2*=(vCoord[0][0]-vCoord[2][0]);
      CT L2=det_inv*(det+det2);
	
      if (L1<0) L1=0; if (L1>1) L1=1;
      if (L2<0) L2=0; if (L2>1) L2=1;
     
      CT L3= CT(1.0) - L1 - L2;
      // Force the point to be inside the simplex if it is not (computation errors)
      if (L3<0)
	{	 	  
	  CT norm_inv=CT(1.0)/(L1+L2);
	  L1*=norm_inv;L2*=norm_inv;L3=0;
	}    
      
      return hlp::numericStaticCast<RT>(L1*vCoord[0][2]+L2*vCoord[1][2]+L3*vCoord[2][2]);
    }

    template <class T, class RT=T, class CT=T>
    static RT interpolateY(const T vCoord[][NDIM], const T pCoord[NDIM])
    {      
      //static const double tol=1.e5*std::numeric_limits<CT>::epsilon();
      static const double tol=1.e4*std::numeric_limits<double>::epsilon();
      CT det=(vCoord[0][0]-vCoord[2][0]);
      det*=(vCoord[1][2]-vCoord[2][2]);
      CT det2=(vCoord[0][2]-vCoord[2][2]);
      det2*=(vCoord[1][0]-vCoord[2][0]);      

      // Check if the determinant is precise enough
      CT det_inv=det-det2;
#if defined(HAVE_BOOST) && defined(HAVE_GMP)      
      if ((fabs(det_inv)<std::max(fabs(det2),fabs(det))*tol)||(det_inv==0))
    	{
	   typedef boost::multiprecision::mpf_float mpfloat;
	   mpfloat a=vCoord[0][0];a-=vCoord[2][0];
	   mpfloat b=vCoord[1][2];b-=vCoord[2][2];	   
	   a*=b;

	   b=vCoord[0][2];b-=vCoord[2][2];	   
	   mpfloat c=vCoord[1][0];c-=vCoord[2][0];	   
	   b*=c;
	   
	   // because of SOS, a degenerate face cannot be intersected, so a!=b
	   c=1;
	   if (a!=b) c/=(a-b); 
	   det_inv = c.convert_to<CT>();	  
	}
      else
	{
	  det=det_inv;
	  det2=1;
	  det_inv=det2/det;	  
	}
#else
      if (det==det2) 
	return hlp::numericStaticCast<RT>((vCoord[0][1]+vCoord[1][1]+vCoord[2][1])/3);
      det=det_inv;
      det2=1;
      det_inv=det2/det;      
#endif

      det=(pCoord[0]-vCoord[2][0]);
      det*=(vCoord[1][2]-vCoord[2][2]);
      det2=(pCoord[2]-vCoord[2][2]);
      det2*=(vCoord[2][0]-vCoord[1][0]);
      CT L1=det_inv*(det+det2);

      det=(pCoord[0]-vCoord[2][0]);
      det*=(vCoord[2][2]-vCoord[0][2]);
      det2=(pCoord[2]-vCoord[2][2]);
      det2*=(vCoord[0][0]-vCoord[2][0]);
      CT L2=det_inv*(det+det2);
      
      if (L1<0) L1=0; if (L1>1) L1=1;
      if (L2<0) L2=0; if (L2>1) L2=1;
     
      CT L3 = CT(1.0) - L1 - L2;
      // Force the point to be inside the simplex if it is not (computation errors)
      if (L3<0)
	{	 	  
	  CT norm_inv=CT(1.0)/(L1+L2);
	  L1*=norm_inv;L2*=norm_inv;L3=0;
	}
     
      return hlp::numericStaticCast<RT>(L1*vCoord[0][1]+L2*vCoord[1][1]+L3*vCoord[2][1]);
    }

    template <class T, class RT=T, class CT=T>
    static RT interpolateX(const T vCoord[][NDIM], const T pCoord[NDIM])
    {
      // int dbg=false;
      // if ((fabs(pCoord[1]-1) < 1.E-8)&&
      // 	  (fabs(pCoord[2]-0.5) < 1.E-8)&&
      // 	  (vCoord[0][0]>0.597)&&(vCoord[0][0]<0.598))
      // 	dbg=true;

      //static const double tol=1.e5*std::numeric_limits<CT>::epsilon();
      static const double tol=1.e4*std::numeric_limits<double>::epsilon();
      CT det=(vCoord[0][1]-vCoord[2][1]);
      det*=(vCoord[1][2]-vCoord[2][2]);
      CT det2=(vCoord[0][2]-vCoord[2][2]);
      det2*=(vCoord[1][1]-vCoord[2][1]);  
    
      // if (dbg) std::cout <<"det=("<<vCoord[0][1]<<"-"<<vCoord[2][1]<<")*("
      // 			 <<vCoord[1][2]<<"-"<<vCoord[2][2]<<")"<<std::endl;
      // if (dbg) std::cout <<"det2=("<<vCoord[0][2]<<"-"<<vCoord[2][2]<<")*("
      // 			 <<vCoord[1][1]<<"-"<<vCoord[2][1]<<")"<<std::endl;
      
      // if (dbg) std::cout <<"det-det2 = "<<det<<"-"<<det2<<"="<<det-det2<<" (<"
      // 			 <<std::max(fabs(det2),fabs(det))*tol<<"?)"<<std::endl;

      // Check if the determinant is precise enough 
      CT det_inv=det-det2;
#if defined(HAVE_BOOST) && defined(HAVE_GMP)            
      if ((fabs(det_inv)<std::max(fabs(det2),fabs(det))*tol)||(det_inv==0))
	{
	  typedef boost::multiprecision::mpf_float mpfloat;
	   mpfloat a=vCoord[0][1];a-=vCoord[2][1];
	   mpfloat b=vCoord[1][2];b-=vCoord[2][2];	   
	   a*=b;

	   b=vCoord[0][2];b-=vCoord[2][2];	   
	   mpfloat c=vCoord[1][1];c-=vCoord[2][1];	   
	   b*=c;
	   
	   // because of SOS, a degenerate face cannot be intersected, so a!=b
	   // FIXME: linear interpolation for truely degenerate cases (in case no SOS is
	   // used ...)
	   c=1;
	   // if (dbg) std::cout << "YES! (a-b)=" << a <<"-"<<b<<"="<<a-b<<std::endl;
	   if (a!=b) c/=(a-b); 
	   det_inv = c.convert_to<CT>();
	}
      else
	{
	  det=det_inv;
	  det2=1;
	  // if (dbg) std::cout << "NO!"<<std::endl;
	  det_inv=det2/det;	  	  
	}      
#else
      //FIXME: linear interpolate ...
      if (det==det2) 
	return hlp::numericStaticCast<RT>((vCoord[0][0]+vCoord[1][0]+vCoord[2][0])/3);
      det=det_inv;
      det2=1;
      det_inv=det2/det;
#endif
      // if (dbg) std::cout << "det_inv = "<<det_inv<<std::endl;
      det=(pCoord[1]-vCoord[2][1]);
      det*=(vCoord[1][2]-vCoord[2][2]);
      det2=(pCoord[2]-vCoord[2][2]);
      det2*=(vCoord[2][1]-vCoord[1][1]);
      CT L1=det_inv*(det+det2);
      // if (dbg) std::cout << "L1 = "<<det_inv<<"*'"<<det<<"+"<<det2<<")="<<L1<<std::endl;
      det  =(pCoord[1]-vCoord[2][1]);
      det *=(vCoord[2][2]-vCoord[0][2]);
      det2 =(pCoord[2]-vCoord[2][2]);
      det2*=(vCoord[0][1]-vCoord[2][1]);
      CT L2=det_inv*(det+det2);
      // if (dbg) std::cout << "L2 = "<<det_inv<<"*'"<<det<<"+"<<det2<<")="<<L2<<std::endl;

      if (L1<0) L1=0; if (L1>1) L1=1;
      if (L2<0) L2=0; if (L2>1) L2=1;
      
      CT L3= CT(1.0) - L1 - L2;
      // if (dbg) std::cout << "L3=" <<  L3<<std::endl;
      // Force the point to be inside the simplex if it is not (computation errors)
      if (L3<0)
	{	 	  
	  CT norm_inv=CT(1.0)/(L1+L2);
	  L1*=norm_inv;L2*=norm_inv;L3=0;
	  // if (dbg) std::cout << "Corrected L3! L1=" 
	  // 		     <<L1<<" L2=" <<L2<<" L3=" <<L3<<std::endl;
	}

      // if (dbg) std::cout << "Result : " 
      // 			 << L1<<"*"<<vCoord[0][0]<<"+"
      // 			 << L2<<"*"<<vCoord[1][0]<<"+"
      // 			 << L3<<"*"<<vCoord[2][0]<<"="
      // 			 <<L1*vCoord[0][0]+L2*vCoord[1][0]+L3*vCoord[2][0]<<"="
      // 			 <<hlp::numericStaticCast<RT>(L1*vCoord[0][0]+L2*vCoord[1][0]+L3*vCoord[2][0])<<std::endl;
      return hlp::numericStaticCast<RT>(L1*vCoord[0][0]+L2*vCoord[1][0]+L3*vCoord[2][0]);
    }
    
    template <class T, class RT=T, class CT=T>
    static RT interpolate(const T vCoord[NDIM][NDIM], const T pCoord[NDIM], int dim)
    {
      if (dim==2) return interpolateZ<T,RT,CT>(vCoord,pCoord);
      else if (dim==1) return interpolateY<T,RT,CT>(vCoord,pCoord);
      else return interpolateX<T,RT,CT>(vCoord,pCoord);
    }
   
    template <class T, class G, class RT=T, class CT=T>
    static RT interpolate(const T vCoord[NDIM][NDIM], const T pCoord[NDIM],
			      int dim, const G* geometry)      
    {
      T checkedVCoord[NDIM][NDIM];
      std::copy_n(&vCoord[0][0],NDIM*NDIM,&checkedVCoord[0][0]);
      geometry->template checkCoordsConsistency<T,NDIM,NDIM>(checkedVCoord,pCoord);

      if (dim==2) return interpolateZ<T,RT,CT>(checkedVCoord,pCoord);
      else if (dim==1) return interpolateY<T,RT,CT>(checkedVCoord,pCoord);
      else return interpolateX<T,RT,CT>(checkedVCoord,pCoord);
    }


  private:
    static bool findVertex(double x1,double y1,double z1,double vx, double vy, double vz, double tol=2.E-5)
    {
      if ((fabs(vx-x1)<=tol)&&
	  (fabs(vy-y1)<=tol)&&
	  (fabs(vz-z1)<=tol))
	{
	  return true;
	}
      return false;
    }

  };

}

#ifdef DEBUG_ME
#undef DEBUG_ME
#endif

#include "../../../internal/namespace.footer"
#endif
