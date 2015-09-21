#ifndef __POINT_IN_SIMPLEX_BASE_3D_HXX__
#define __POINT_IN_SIMPLEX_BASE_3D_HXX__

#include <string.h> // for memcpy

#include "pointInSimplexBasePrototype.hxx"
#include "pointInSimplexBase_2D.hxx"

#include "../predicates/orientation.hxx"

#include "../../geometry/simplexInterpolator.hxx"
#include "../../geometry/geometricProperties.hxx"

#include "../../tools/wrappers/boostMultiprecisionFloat128.hxx"

#ifdef HAVE_BOOST
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#endif
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#endif

#include "../../internal/namespace.header"

//#define DEBUG_ME 1

namespace internal {  
  
  template <int filterType>
  class PointInSimplexBaseT<3,filterType>
  {

  public:
    static const int NDIM=3;
    typedef predicate::OrientationT<3,filterType> Orientation;
   
    template <class T, class T2, class T3>
    static int getRayIntersection_OLD(const T vCoord[][NDIM], const T pCoord[NDIM], 
				  const int dim, T2 faceIndex[2], T3 zCoord[2])
    {
#ifdef DEBUG_ME
      bool debug=false;
      const double findCoord[3]={-0.101562,0.101562,0.101562};
      if (findVertex(pCoord[0],pCoord[1],pCoord[2],findCoord[0],findCoord[1],findCoord[2])) 
	{
	  debug=true;
	  printf("Debugging point (%g %g %g) with simplex:\n",pCoord[0],pCoord[1],pCoord[2]);
	  for (int i=0;i<NDIM+1;++i) 
	    printf("    P[%d]=(%g %g %g)\n",i,
		   vCoord[i][0],vCoord[i][1],vCoord[i][2]);

	}
#endif

      int coord[2];
      T vCoord2D[NDIM][2];
      T pCoord2D[2];
      T faceCoord[2][NDIM][NDIM];     
      int nFaces=0;
      
      int j=0;
      for (int i=0;i<NDIM;++i) 
	{
	  if (i!=dim) 
	    {
	      coord[j]=i;
	      pCoord2D[j]=pCoord[i];
	      ++j;
	    }	 
	}

      // Now count, when projected along Z axis, how many of the 4 triangles contain P.
      // Answer may only be 2 or 0 !
      // NOTE: the order must correspond to index2D[4][3]
      vCoord2D[0][0]=vCoord[0][coord[0]];vCoord2D[0][1]=vCoord[0][coord[1]];
      vCoord2D[1][0]=vCoord[1][coord[0]];vCoord2D[1][1]=vCoord[1][coord[1]];
      vCoord2D[2][0]=vCoord[2][coord[0]];vCoord2D[2][1]=vCoord[2][coord[1]];
      if (test2D<T>(vCoord2D,pCoord2D))
	{
	  memcpy(faceCoord[nFaces][0],vCoord[0],sizeof(T)*NDIM*3);
	  faceIndex[nFaces]=3;
	  nFaces++;
#ifdef DEBUG_ME
	  if (debug) printf("FOUND (%g %g) @3 : (%g %g) (%g %g) (%g %g)\n",
			    pCoord2D[0],pCoord2D[1],vCoord2D[0][0],vCoord2D[0][1],
			    vCoord2D[1][0],vCoord2D[1][1],vCoord2D[2][0],vCoord2D[2][1]);
#endif
	}

      vCoord2D[0][0]=vCoord[1][coord[0]];vCoord2D[0][1]=vCoord[1][coord[1]];
      vCoord2D[1][0]=vCoord[2][coord[0]];vCoord2D[1][1]=vCoord[2][coord[1]];
      vCoord2D[2][0]=vCoord[3][coord[0]];vCoord2D[2][1]=vCoord[3][coord[1]];
      if (test2D<T>(vCoord2D,pCoord2D))
	{
	  memcpy(faceCoord[nFaces][0],vCoord[1],sizeof(T)*NDIM*3);
	  faceIndex[nFaces]=0;
	  nFaces++;
#ifdef DEBUG_ME
	  if (debug) printf("FOUND (%g %g) @0 : (%g %g) (%g %g) (%g %g)\n",
			    pCoord2D[0],pCoord2D[1],vCoord2D[0][0],vCoord2D[0][1],
			    vCoord2D[1][0],vCoord2D[1][1],vCoord2D[2][0],vCoord2D[2][1]);
#endif
	}
      
      // At most 2 faces may be crossed
      if (nFaces<2)
	{
	  vCoord2D[0][0]=vCoord[0][coord[0]];vCoord2D[0][1]=vCoord[0][coord[1]];
	  vCoord2D[1][0]=vCoord[2][coord[0]];vCoord2D[1][1]=vCoord[2][coord[1]];
	  vCoord2D[2][0]=vCoord[3][coord[0]];vCoord2D[2][1]=vCoord[3][coord[1]];
	  if (test2D<T>(vCoord2D,pCoord2D))
	    {
	      memcpy(faceCoord[nFaces][0],vCoord[0],sizeof(T)*NDIM);
	      memcpy(faceCoord[nFaces][1],vCoord[2],sizeof(T)*NDIM*2);
	      faceIndex[nFaces]=1;
	      nFaces++;
#ifdef DEBUG_ME
	      if (debug) printf("FOUND (%g %g) @1 : (%g %g) (%g %g) (%g %g)\n",
				pCoord2D[0],pCoord2D[1],vCoord2D[0][0],vCoord2D[0][1],
				vCoord2D[1][0],vCoord2D[1][1],vCoord2D[2][0],vCoord2D[2][1]);
#endif
	    }	  

	  // at least 3 faces are not crossed => 0 are ...
	  if (nFaces==1)
	    {
	      vCoord2D[0][0]=vCoord[0][coord[0]];vCoord2D[0][1]=vCoord[0][coord[1]];
	      vCoord2D[1][0]=vCoord[1][coord[0]];vCoord2D[1][1]=vCoord[1][coord[1]];
	      vCoord2D[2][0]=vCoord[3][coord[0]];vCoord2D[2][1]=vCoord[3][coord[1]];
	      if (test2D<T>(vCoord2D,pCoord2D))
		{
		  memcpy(faceCoord[nFaces][0],vCoord[0],sizeof(T)*NDIM);
		  memcpy(faceCoord[nFaces][1],vCoord[1],sizeof(T)*NDIM);
		  memcpy(faceCoord[nFaces][2],vCoord[3],sizeof(T)*NDIM);
		  faceIndex[nFaces]=2;
		  nFaces++;
#ifdef DEBUG_ME
		  if (debug) printf("FOUND (%g %g) @2 : (%g %g) (%g %g) (%g %g)\n",
				    pCoord2D[0],pCoord2D[1],vCoord2D[0][0],vCoord2D[0][1],
				    vCoord2D[1][0],vCoord2D[1][1],vCoord2D[2][0],vCoord2D[2][1]);
#endif
		}
	    }
	}

      if (nFaces==2)
	{	  
	  // Sort the two faces so that the points are always in the same order, whatever 
	  // the original order of the simplices. This ensure IEEE FP errors consistency !
	  for (int i=0;i<2;++i)
	    for (int j=0;j<NDIM-1;++j)
	      for (int k=j+1;k<NDIM;++k)
		{
		  if (faceCoord[i][j][0]<faceCoord[i][k][0])
		    for (int l=0;l<NDIM;++l) std::swap(faceCoord[i][j][l],faceCoord[i][k][l]);
		  else if (faceCoord[i][j][0]==faceCoord[i][k][0])
		    {
		      if (faceCoord[i][j][1]<faceCoord[i][k][1])
			for (int l=0;l<NDIM;++l) std::swap(faceCoord[i][j][l],faceCoord[i][k][l]);
		      else if (faceCoord[i][j][1]==faceCoord[i][k][1])
			{
			  if (faceCoord[i][j][2]<faceCoord[i][k][2])
			    for (int l=0;l<NDIM;++l) std::swap(faceCoord[i][j][l],faceCoord[i][k][l]);
			}
		    }
		}
	  
	  if (dim==2)
	    {
	      zCoord[0]=interpolateZ<T>(faceCoord[0],pCoord);
	      zCoord[1]=interpolateZ<T>(faceCoord[1],pCoord);		  
	    }
	  else if (dim==1)
	    {
	      zCoord[0]=interpolateY<T>(faceCoord[0],pCoord);
	      zCoord[1]=interpolateY<T>(faceCoord[1],pCoord);
	    }
	  else
	    {
	      zCoord[0]=interpolateX<T>(faceCoord[0],pCoord);
	      zCoord[1]=interpolateX<T>(faceCoord[1],pCoord);
	    }
	  
#ifdef DEBUG_ME
	  if (debug) printf("limit@%d : %.17g < %.17g < %.17g (%d)\n",dim,
			    zCoord[0],pCoord[dim],zCoord[1],
			    ((zCoord[0]<pCoord[dim])!=(zCoord[1]<pCoord[dim])));
#endif

	  return 2;
	  //return ((z0<pCoord[2])!=(z1<pCoord[2]));	  
	}
      
      return 0;
    }

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

    // FIXME: IT IS UGLY TO DUPLICATE THIS FUNCTION, FIND A BETTER WAY ?
    template <class T, class T2, class CT=T2>
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
	      intersectionCoord = interpolate<T,T2,CT>(vCoord,pCoord,dim);		
	      return 1;
	    }
	}
	
      return 0;
    }

    template <class T, class T2, class G, class CT=T2>
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
      
      /*
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
      */

      /*DELETEME*/
      /*
      bool test=false;
      double vertC[2][NDIM]={{-0.015625,0.015625,-0.273438},
			    {-0.015625,0.015625,-0.265625}};
      if (findVertex(pCoord[0],pCoord[1],pCoord[2],
		     vertC[0][0],vertC[0][1],vertC[0][2])||
	  findVertex(pCoord[0],pCoord[1],pCoord[2],
		     vertC[1][0],vertC[1][1],vertC[1][2]))
	test=true;
      if (test)
	printf("Testing point @(%20.20e %20.20e) along d%d\n   -> in (%20.20e %20.20e)\n         (%20.20e %20.20e)\n         (%20.20e %20.20e)\n",
	       pCoord2D[0],pCoord2D[1],dim,	     
	       vCoord2D[0][0],vCoord2D[0][1],
	       vCoord2D[1][0],vCoord2D[1][1],
	       vCoord2D[2][0],vCoord2D[2][1]);
	       */

      // FIXME: This will only work for a CUBIC box as geometry is 3D => coords always correspond to the first 2 dims !
      if (test2D<T>(vCoord2D,pCoord2D,dim,geometry))
	  {
	    double testPointCoord[NDIM];
	    std::copy(pCoord,pCoord+NDIM,testPointCoord);
	    int res1 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
					 testPointCoord,geometry);
      
	    testPointCoord[dim]=otherCoord;
	    int res2 = Orientation::test(vCoord[0],vCoord[1],vCoord[2],
					 testPointCoord,geometry);
	    /*DELETEME*/ //if (test) printf("Passed test, (%d/%d)\n",res1,res2);
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
		intersectionCoord = interpolate<T,G,T2,CT>(vCoord,pCoord,dim,geometry);
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
      
      /*
      int faceIndex[2];
      double dimCoord[2];
      static const int dim=2;
      int result = getRayIntersection(vCoord, pCoord,dim,faceIndex,dimCoord);
				      
      if (result==2)
	return ((dimCoord[0]<pCoord[dim])!=(dimCoord[1]<pCoord[dim]));
      return 0;
      */
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
    /*
    template <class T>
    static int test2D_X(const T vCoord[][2], const T pCoord[NDIM-1])
    {
      return PointInSimplexBaseT<2,intersectionType>::template testX<T>(vCoord,pCoord);
    }

     template <class T>
    static int test2D_Y(const T vCoord[][2], const T pCoord[NDIM-1])
    {
      return PointInSimplexBaseT<2,intersectionType>::template testY<T>(vCoord,pCoord);
    }
    */

    /*
    template <class T, class T2>
    static int planeSide(const T vCoord[][NDIM], const T* const pCoord,const T2* const index)
    {    
      double det1 = 
	(vCoord[1][1]-vCoord[0][1])*(vCoord[2][2]-vCoord[0][2])-
	(vCoord[2][1]-vCoord[0][1])*(vCoord[1][2]-vCoord[0][2]);
      double det2 = 
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][2]-vCoord[0][2])-
	(vCoord[2][0]-vCoord[0][0])*(vCoord[1][2]-vCoord[0][2]);
      double det3 = 
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][1]-vCoord[0][1])-
	(vCoord[2][0]-vCoord[0][0])*(vCoord[1][1]-vCoord[0][1]);

      double det = det1*(pCoord[0]-vCoord[0][0])-det2*(pCoord[1]-vCoord[0][1])+det3*(pCoord[2]-vCoord[0][2]);
      
      // det is NULL => simulation of simplicity
      if (det==0)
	{
	  det=-det3;
	  if (det==0) 
	    {
	      det=det2; // This is indeed a + sign (sign of epsilon times -1)
	      if (det==0) 
		{
		  det=-det1;
		  if (det==0)
		    {
		      glb::console->print<LOG_WARNING>("SIMULATING SIMPLICITY\n");
		      det1 = 
			(pCoord[1]-vCoord[0][1])*(vCoord[1][0]-vCoord[2][0])-
			(pCoord[0]-vCoord[0][0])*(vCoord[1][1]-vCoord[0][1]);
		      det2 =
			(pCoord[1]-vCoord[0][1])*(vCoord[2][0]-vCoord[0][0])-
			(pCoord[0]-vCoord[0][0])*(vCoord[2][1]-vCoord[0][1]);
		      det3 =
			(pCoord[0]-vCoord[0][0])*(vCoord[1][1]-vCoord[0][1])-
			(pCoord[1]-vCoord[0][1])*(vCoord[1][0]-vCoord[0][0]);
		      
		      if ((index[0]<index[1])&&(index[0]<index[2])) 
			{det=det1;det1=det2;det2=det3;index[0]=index[1];index[1]=index[2];}
		      else if (index[1]<index[2]) 
			{det=det2;det2=det3;index[1]=index[2];}
		      else 
			{det=det3;}
		      if (det==0)
			{
			  if (index[0]<index[1]) 
			    {det=det1;det1=det2;index[0]=index[1];}
			  else 
			    {det=det2;}
			  if (det==0) 
			    {
			      det=det1;
			    }
			}
		    }
		}
	    }
	}
      if (det==0) glb::console->print<LOG_WARNING>("NULL DETERMINANT, SOS failed \n");
      return (det>0);     
    }
    */

    /*
    template <class T>
    static int planeSide(const T vCoord[][NDIM], const T pCoord[NDIM])
    {
      double det1 = 
	(vCoord[1][1]-vCoord[0][1])*(vCoord[2][2]-vCoord[0][2])-
	(vCoord[2][1]-vCoord[0][1])*(vCoord[1][2]-vCoord[0][2]);
      double det2 = 
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][2]-vCoord[0][2])-
	(vCoord[2][0]-vCoord[0][0])*(vCoord[1][2]-vCoord[0][2]);
      double det3 = 
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][1]-vCoord[0][1])-
	(vCoord[2][0]-vCoord[0][0])*(vCoord[1][1]-vCoord[0][1]);

      double det = det1*(pCoord[0]-vCoord[0][0])-det2*(pCoord[1]-vCoord[0][1])+det3*(pCoord[2]-vCoord[0][2]);
      
      // Simulation of simplicity (only pCoord is perturbed, or equivalently works if vCoord perturbation is smaller)
      if (det==0)
	{
	  det=-det3;
	  if (det==0) 
	    {
	      det=det2; // This is indeed a + sign (sign of epsilon times -1)
	      if (det==0) 
		{
		  det=-det1;
		}
	    }
	  if (det==0) glb::console->print<LOG_WARNING>("NULL DETERMINANT, SOS failed \n");
	}
      
      return (det>0);   
    }
    */

    /*
    template <class T>
    static double interpolateZ(const T vCoord[NDIM][NDIM], const T pCoord[NDIM-1])
    {
      const T z[NDIM]={vCoord[0][2],vCoord[1][2],vCoord[2][2]};
      const T v[NDIM][NDIM-1]={vCoord[0][0],vCoord[0][1],vCoord[1][0],vCoord[1][1],vCoord[2][0],vCoord[2][1]};
      SimplexInterpolatorT<2> sI(vCoord);
      return sI.template interpolateInside<T,T,1>(pCoord,z);
    }
    */
    
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

    //#ifdef DEBUG_ME
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
    //#endif
  };

}

#ifdef DEBUG_ME
#undef DEBUG_ME
#endif

#include "../../internal/namespace.footer"
#endif
