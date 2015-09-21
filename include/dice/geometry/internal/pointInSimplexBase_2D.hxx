#ifndef __POINT_IN_SIMPLEX_BASE_2D_HXX__
#define __POINT_IN_SIMPLEX_BASE_2D_HXX__

#include "pointInSimplexBasePrototype.hxx"
#include "../predicates/orientation.hxx"

#include "../../internal/namespace.header"

namespace internal {
  
  template <int filterType>
  class PointInSimplexBaseT<2,filterType>
  {
  public:
    static const int NDIM=2;
    typedef predicate::OrientationT<2,filterType> Orientation;
    /*
    template <class T, class T2, class T3>
    static int getRayIntersection(const T vCoord[][NDIM], const T pCoord[NDIM], 
				  const int dim, T2 faceIndex[2], T3 zCoord[2])
    {
      printf("SORRY: getRayIntersection is not implemented in 2D ! (but that's easy ;) )\n");
      exit(-1);
    }
    */
    template <class T, class T2>
    static int getRayFacetIntersection(T vCoord[][NDIM], const T pCoord[NDIM], 
				       const int dim, T2 &zCoord)
    {
      printf("SORRY: getRayIntersection is not implemented in 2D ! (but that's easy ;) )\n");
      exit(-1);
    }

    template <class T, class T2, class G>
    static int getRayFacetIntersection(T vCoord[][NDIM], const T pCoord[NDIM], 
				       const int dim, T2 &zCoord, const G* geometry)
    {
      printf("SORRY: getRayIntersection is not implemented in 2D ! (but that's easy ;) )\n");
      exit(-1);
    }

    template <class T, class T2, class CT=T2>
    static int getSegmentFacetIntersection(const T vCoord[][NDIM], const T pCoord[NDIM], 
					   int dim, T p2, T2 &zCoord)
    {
      printf("SORRY: getRayIntersection is not implemented in 2D ! (but that's easy ;) )\n");
      exit(-1);
    }

    template <class T, class T2, class G, class CT=T2>
    static int getSegmentFacetIntersection(const T vCoord[][NDIM], const T pCoord[NDIM], 
					   int dim, T p2, T2 &zCoord, const G* geometry)
    {
      printf("SORRY: getRayIntersection is not implemented in 2D ! (but that's easy ;) )\n");
      exit(-1);
    }
   

    template <class T>
    static int test(const T vCoord[][NDIM], const T* const pCoord)
    {
      //return testY(vCoord,pCoord);
      
      int res=0;
      res+=Orientation::test(vCoord[0],vCoord[1],pCoord);
      res+=Orientation::test(vCoord[1],vCoord[2],pCoord);
      res+=Orientation::test(vCoord[2],vCoord[0],pCoord);
      return ((res==3)||(res==0));      
    }

    template <class T, class G>
    static int test(const T vCoord[][NDIM], const T* const pCoord, const G* geometry)
    {
      //return testY(vCoord,pCoord);
      
      int res=0;
      res+=Orientation::test(vCoord[0],vCoord[1],pCoord,geometry);
      res+=Orientation::test(vCoord[1],vCoord[2],pCoord,geometry);
      res+=Orientation::test(vCoord[2],vCoord[0],pCoord,geometry);
      /*
      if ((fabs(pCoord[0]-(-6.2890625E-1))<1.E-5)&&
	  (fabs(pCoord[1]-(1))<1.E-5))
	{
	  printf("TESTING: (%20.20e %20.20e) (%20.20e %20.20e) (%20.20e %20.20e) -> (%d,%d,%d)=>%d\n",
		 vCoord[0][0],vCoord[0][1],vCoord[1][0],vCoord[1][1],vCoord[2][0],vCoord[2][1],
		 Orientation::test(vCoord[0],vCoord[1],pCoord,geometry),
		 Orientation::test(vCoord[1],vCoord[2],pCoord,geometry),
		 Orientation::test(vCoord[2],vCoord[0],pCoord,geometry),
		 ((res==3)||(res==0)));
	  
	}
      */
      return ((res==3)||(res==0));      
    }

  private:
    
    // Check whether a point is inside or outside a simplex with vertices vCoord. This also
    // work for any polygon with any number of vertices.
    // This will consistently return true or false, no boundary here ...
    // This is adapted from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
    // with some corrections to make it robust wrt IEEE rounding errors.
    template <class T>
    static int testX(const T vCoord[][NDIM], const T* const pCoord)
    {
      int i, j, c = 0;
      for (i = 0, j = NDIM; i < NDIM+1; j = i++) {
	if  ((vCoord[i][1]>pCoord[1]) != (vCoord[j][1]>pCoord[1]))
	  {	  
	    if (vCoord[j][0]<vCoord[i][0])
	      {
		if ( pCoord[0] < 
		     (vCoord[j][0]-vCoord[i][0])/(vCoord[j][1]-vCoord[i][1])*
		     (pCoord[1]-vCoord[i][1])  + vCoord[i][0] )
		  {c=!c;}
	      }
	    else if ( pCoord[0] < 
		      (vCoord[i][0]-vCoord[j][0])/(vCoord[i][1]-vCoord[j][1])*
		      (pCoord[1]-vCoord[j][1])  + vCoord[j][0] )
	      {c=!c;}
	  }	    
      }
      return c;
    }
    
    // Check whether a point is inside or outside a simplex with vertices vCoord. This also
    // work for any polygon with any number of vertices.
    // This will consistently return true or false, no boundary here ...
    // This is adapted from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
    // with some corrections to make it robust wrt IEEE rounding errors, with ray 
    // along Y instead of X, and with a different perturbation for SoS to comply to our 
    // implementation.
    template <class T>
    static int testY(const T vCoord[][NDIM], const T* const pCoord)
    {
      int i, j, c = 0;
      for (i = 0, j = NDIM; i < NDIM+1; j = i++) {
	if  ((vCoord[i][0]<pCoord[0]) != (vCoord[j][0]<pCoord[0]))
	  {	  
	    if (vCoord[j][1]<vCoord[i][1])
	      {
		if ( pCoord[1] <=
		     (vCoord[j][1]-vCoord[i][1])/(vCoord[j][0]-vCoord[i][0])*
		     (pCoord[0]-vCoord[i][0])  + vCoord[i][1] )
		  {c=!c;}
	      }
	    else if ( pCoord[1] <=
		      (vCoord[i][1]-vCoord[j][1])/(vCoord[i][0]-vCoord[j][0])*
		      (pCoord[0]-vCoord[j][0])  + vCoord[j][1] )
	      {c=!c;}
	  }	    
      }
      return c;
    }

    /*
public:
    // Check whether a point is inside or outside a simplex with vertices vCoord. This also
    // work for any polygon with any number of vertices.
    // This will consistently return true or false, no boundary here ...
    // from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
    template <class T>
    static int test11(const T vCoord[][NDIM], const T* const pCoord)
    {
      int i, j, c = 0;
      for (i = 0, j = NDIM; i < NDIM+1; j = i++) {
	if  ((vCoord[i][1]>pCoord[1]) != (vCoord[j][1]>pCoord[1]))
	  {
	    if (pCoord[1]!=vCoord[j][1])
	      {
		
		// if (((vCoord[j][0]-vCoord[i][0])/(vCoord[j][1]-vCoord[i][1])*
		//      (pCoord[1]-vCoord[i][1])  + vCoord[i][0]) != 
		//     ((vCoord[j][0]-vCoord[i][0])*(pCoord[1]-vCoord[i][1])/
		//      (vCoord[j][1]-vCoord[i][1]) + vCoord[i][0]))
		       
		//   {
		//     printf("FIDD\n");
		//     printf("%.18lg/%.18lg*%.18lg+%.18lg = %.18lg != %.18lg\n",
		// 	   (vCoord[j][0]-vCoord[i][0]),
		// 	   (vCoord[j][1]-vCoord[i][1]),
		// 	   (pCoord[1]-vCoord[i][1]),
		// 	   vCoord[i][0],
		// 	   (vCoord[j][0]-vCoord[i][0])/(vCoord[j][1]-vCoord[i][1])*
		// 	   (pCoord[1]-vCoord[i][1])  + vCoord[i][0],
		// 	   (vCoord[j][0]-vCoord[i][0])*(pCoord[1]-vCoord[i][1])/
		// 	   (vCoord[j][1]-vCoord[i][1]) + vCoord[i][0]);
		//     exit(-1);
		//   }
		
	      }

	    if (pCoord[1]==vCoord[j][1])
	      {
		if (pCoord[0]<vCoord[j][0]) c = !c;
	      }
	    else if ( pCoord[0] < 
		      (vCoord[j][0]-vCoord[i][0])/(vCoord[j][1]-vCoord[i][1])*
		      (pCoord[1]-vCoord[i][1])  + vCoord[i][0] )
	      {
		c = !c;
	      }
	  }
	    
      }
      return c;
    }
    
    // Check wether a point is inside or outside a simplex with vertices vCoord. This also
    // work for any polygon with any number of vertices.
    // This will consistently return true or false, no boundary here ...
    // from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
    template <class T>
    static int test22(const T vCoord[][NDIM], const T* const pCoord)
    {
      int i, j, c = 0;
      for (i = 0, j = NDIM; i < NDIM+1; j = i++) {
	if ( ((vCoord[i][1]>pCoord[1]) != (vCoord[j][1]>pCoord[1])) &&
	     ( 
	      ((pCoord[1]==vCoord[j][1])&&(pCoord[0]<vCoord[j][0])) ||
	      (pCoord[0] < (vCoord[j][0]-vCoord[i][0])/(vCoord[j][1]-vCoord[i][1])*(pCoord[1]-vCoord[i][1])  + vCoord[i][0])
	       ) 
	     )
	  c = !c;
      }
      return c;
    }

    static int pnpolyC(int nvert, double *vertx, double *verty, double testx, double testy)
    {
      int i, j, c = 0;
      for (i = 0, j = nvert-1; i < nvert; j = i++) {
	if ((verty[i]>testy) != (verty[j]>testy))
	  {
	    printf("Test1[%d,%d]= (%d != %d) == 1\n",i,j,(verty[i]>testy),(verty[j]>testy));
	    printf("(testy==verty[%d]) == (%.18lg==%.18lg) = %d.\n",
		   j,testy,verty[j],(testy==verty[j]));
	    if (testy!=verty[j])
	      printf("Testing %.18lg<(%.18lg-%.18lg)/(%.18lg-%.18lg)*(%.18lg-%.18lg)+%.18lg\n --> %.18lg<%.18lg/%.18lg*%.18lg\n",
		     testx,vertx[j],vertx[i],verty[j],verty[i],testy,verty[i],vertx[i],
		     testx-vertx[i],vertx[j]-vertx[i],verty[j]-verty[i],testy-verty[i]);
	    
	    if (testy==verty[j])
	      {
		printf("(testx<vertx[%d]) == (%.18lg<%.18lg) = %d.\n",
		       j,testx,vertx[j],(testx<vertx[j]));
		if (testx<vertx[j])
		  {
		    printf("-----> c=%d => %d\n",c,!c);
		    c=!c;
		  }
		else printf("-----> c=%d\n",c);
	      }
	    else if (testx < (vertx[j]-vertx[i])/ (verty[j]-verty[i]) * (testy-verty[i])  + vertx[i])
	      {
		printf("-----> c=%d => %d\n",c,!c);
		c=!c;
	      }
	    else printf("-----> c=%d\n",c);
	  }	
	else printf("Test1[%d,%d]= (%d != %d) == 0\n",i,j,(verty[i]>testy),(verty[j]>testy));
	  
      }
      return c;
    }

    static int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
    {
      int i, j, c = 0;
      for (i = 0, j = nvert-1; i < nvert; j = i++) {
	if ( ((verty[i]>testy) != (verty[j]>testy)) )
	  {
	    printf("Test1[%d,%d]= (%d != %d) == 1\n",i,j,(verty[i]>testy),(verty[j]>testy));
	    printf("Testing %.18lg<(%.18lg-%.18lg)*(%.18lg-%.18lg)/(%.18lg-%.18lg)+%.18lg\n --> %.18lg<(%.18lg*%.18lg)/%.18lg\n",
		   testx,vertx[j],vertx[i],testy,verty[i],verty[j],verty[i],vertx[i],
		   testx-vertx[i],vertx[j]-vertx[i],testy-verty[i],verty[j]-verty[i]);
	    if ((testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]))
	      printf("-----> c=%d => %d\n",c,!c);
	    else
	      printf("-----> c=%d\n",c);
	  }
	else
	  printf("Test1[%d,%d]= (%d != %d) == 0\n",i,j,(verty[i]>testy),(verty[j]>testy));


	if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
	  c = !c;
      }
      return c;
    }
*/

  };

} // namespace internal
 
#include "../../internal/namespace.footer"
#endif
