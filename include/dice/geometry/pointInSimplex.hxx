#ifndef __POINT_IN_SIMPLEX_HXX__
#define __POINT_IN_SIMPLEX_HXX__

#include "../geometry/predicates/filterType.hxx"
#include "./internal/pointInSimplexBase_2D.hxx"
#include "./internal/pointInSimplexBase_3D.hxx"

/**
 * @file 
 * @brief  Definition of a class to test whether a point is inside our outside a simplex
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Geometry
 *   \{
 */

/**
 * \class PointInSimplexT
 * \brief A class to test whether a point is inside our outside a simplex or a poly.
 * The test is robust in the sense that if a point was tested against all simplices 
 * in simplicial tesselation of space, it would only test positive with one simplex, even 
 * in the degenerate case where it falls precisely on a simplex boundary. This is achieved 
 * by simulation of simplicity, with the convention that the test point coordinate Pi 
 * along dimension i (0<=i<NDIM) is disturbed as follows : 
 *              Pi -> Pi - pow(eps,NDIM-i)
 *
 * \warning testing whether the queried point is within the bounding box of the simplex
 *  before querying is the responsibility of the user, so not doing it may result in a waste
 *  of computational time ...
 * 
 * \tparam ND The number of dimensions
 * \tparam filterType The type of predicate to use, see predicate::filterType::Type 
 *  Default value is predicate::filterType::Adaptive
 * \todo Embed in a higher dimensional space 
 */

template <int ND, int filterType=predicate::filterType::Adaptive>
class PointInSimplexT : protected internal::PointInSimplexBaseT<ND,filterType>
{
protected:
  typedef internal::PointInSimplexBaseT<ND,filterType> Base;
public:
  static const int NDIM = ND;

  /**
   * \brief Test if a point is inside or outside a simplex. Note that the point cannot be 
   *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
   * \param simplex A pointer to a simplex
   * \param pCoord A pointer to the ND coordinates of the point to test
   * \param geometry A pointer to the geometry of the mesh
   * \param coordsAreConsistent For periodic boundaries, set to 1 if coordinates are
   * consistent (i.e. on the same side of the periodic box), 0 if they are not and -1 if 
   * you do not know. This is useless for non periodic boundaries.
   * \tparam S The Simplex class
   * \tparam G The class describing geometric properties of the bounding box. See e.g.
   *  GeometricPropertiesT
   * \return true if the point with coords \a pCoord is INSIDE the simplex, false otherwise
   */
  template <class S, class G>
  static int test(const S *simplex, const typename S::Coord pCoord[NDIM], 
		  const G* geometry, int coordsAreConsistent=-1)
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT][NDIM];
    for (int n=0;n<S::NVERT;++n)
      for (int d=0;d<NDIM;++d)
	vCoord[n][d]=simplex->getVertex(n)->getCoord(d);
    
    //geometry->template checkCoordsConsistency<Coord,S::NVERT,NDIM>(vCoord,pCoord);
    if ((coordsAreConsistent>0) ||	
	((coordsAreConsistent<0)&&
	 geometry->template coordsAreConsistent<Coord,S::NVERT,NDIM>(vCoord,pCoord))
	)
      return Base::test(vCoord,pCoord);
    else
      return Base::test(vCoord,pCoord,geometry);
    /*
    double vx[3];
    double vy[3];
    double px;
    double py;
    for (int i=0;i<3;++i) vx[i]=simplex->getVertex(i)->getCoord(0);
    for (int i=0;i<3;++i) vy[i]=simplex->getVertex(i)->getCoord(1);
    px=pCoord[0];py=pCoord[1];
    return Base::pnpoly(3,vx,vy,px,py);
    */
    
  }
  
 
  /**
   * \brief Test if a point is inside or outside a simplex. Note that the point cannot be 
   *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
   * \param vCoord (NDIM+1) pointers the (NDIM) coordinates of each vertex of the simplex
   * \param pCoord A pointer to the ND coordinates of the point to test
   * \tparam T The type of the coordinates
   * \return True if the point with coords \a pCoord is INSIDE the simplex, false otherwise
   */
  template <class T>
  static int test(const T (&vCoord)[NDIM+1][NDIM], const T pCoord[NDIM])
  {
    return Base::test(vCoord,pCoord);
  }

  /**
   * \brief Test if a point is inside or outside a simplex. Note that the point cannot be 
   *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
   * \param vCoord the (NDIM) coodinates of the (NDIM+1) vertices of the simplex, and the 
   * coordinates of the point to test in vCoord[NDIM+1][]
   * \tparam T The type of the coordinates
   * \return True if the point with coords \a pCoord is INSIDE the simplex, false otherwise
   */
  template <class T>
  static int test(const T (&vCoord)[NDIM+2][NDIM])
  {
    return Base::test(vCoord,&vCoord[NDIM+1][0]);
  }
  /*
  template <class S, class G, class T2, class T3>
  static int intersectRay(const S *simplex, const typename S::Coord pCoord[NDIM], 
			  const G* geometry,int dim, T2 faceIndex[2], T3 crossCoord[2])
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT][NDIM];
    for (int n=0;n<S::NVERT;++n)
      for (int d=0;d<NDIM;++d)
	vCoord[n][d]=simplex->getVertex(n)->getCoord(d);
    
    geometry->template checkCoordsConsistency<Coord,S::NVERT,NDIM>(vCoord,pCoord);
      
    return Base::getRayIntersection(vCoord,pCoord,dim,faceIndex,crossCoord);
  }

  template <class S, class T2, class T3>
  static int intersectRay(const S *simplex, const typename S::Coord pCoord[NDIM], 
			  int dim, T2 faceIndex[2], T3 crossCoord[2])
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT][NDIM];
    for (int n=0;n<S::NVERT;++n)
      for (int d=0;d<NDIM;++d)
	vCoord[n][d]=simplex->getVertex(n)->getCoord(d);
        
    return Base::getRayIntersection(vCoord,pCoord,dim,faceIndex,crossCoord);
  }
  */

  /*
  template <class S, class G, class T>
  static int intersectRayFacet(const S *simplex, const typename S::Coord pCoord[NDIM], 
			       const G* geometry, int dim, int index, T &intersectionCoord)
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT-1][NDIM];
    int i=0;
    for (int n=0;n<S::NVERT;++n)
      {
	if (n!=index)
	  {
	    for (int d=0;d<NDIM;++d)
	      vCoord[i][d]=simplex->getVertex(n)->getCoord(d);
	    ++i;
	  }
      }
    
    //geometry->template checkCoordsConsistency<Coord,S::NVERT-1,NDIM>(vCoord,pCoord);
    if (geometry->template coordsAreConsistent<Coord,S::NVERT-1,NDIM>(vCoord,pCoord))
      return Base::getRayFacetIntersection(vCoord,pCoord,dim,intersectionCoord);  
    else
      return Base::getRayFacetIntersection(vCoord,pCoord,dim,intersectionCoord,geometry);  
  }

  template <class S, class G, class T>
  static int intersectRayFacet(const S *simplex, const typename S::Coord pCoord[NDIM],
			       int dim, int index, T &intersectionCoord)
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT-1][NDIM];
    int i=0;
    for (int n=0;n<S::NVERT;++n)
      {
	if (n!=index)
	  {
	    for (int d=0;d<NDIM;++d)
	      vCoord[i][d]=simplex->getVertex(n)->getCoord(d);
	    ++i;
	  }
      }
    
    return Base::getRayFacetIntersection(vCoord,pCoord,dim,intersectionCoord);
  }
  */

  template <class F, class G, class T, class CT=T>
  static int intersectSegmentFacet(const F &facet, 
				   const typename F::Coord pCoord[NDIM], 
				   int dim, typename F::Coord otherCoord,
				   const G* geometry, T &intersectionCoord,
				   int coordsAreConsistent=-1)
  {
    /*
    typedef typename F::Coord  Coord;
    typedef typename F::Vertex Vertex;

    Vertex *vertex[F::NVERT];
    Coord vCoord[F::NVERT][NDIM];

    facet.getVertices(vertex);        
    for (int n=0;n<F::NVERT;++n)
      vertex[n]->getCoords(vCoord[n]);
*/
    typedef typename F::Coord Coord;

    Coord vCoord[F::NVERT][NDIM];
    for (int n=0;n<F::NVERT;++n)
      {
	typename F::Vertex *v=facet.getVertex(n);
	for (int d=0;d<NDIM;++d)
	  vCoord[n][d]=v->getCoord(d);
      }
    
    //geometry->template checkCoordsConsistency<Coord,F::NVERT,NDIM>(vCoord,pCoord);
    if ((coordsAreConsistent>0)|| 
	((coordsAreConsistent<0)&& 
	 geometry->template coordsAreConsistent<Coord,F::NVERT,NDIM>(vCoord,pCoord))
	)
      return Base::template getSegmentFacetIntersection<Coord,T,CT>
	(vCoord,pCoord,dim,otherCoord,intersectionCoord);
    else
      return Base::template getSegmentFacetIntersection<Coord,T,G,CT>
	(vCoord,pCoord,dim,otherCoord,intersectionCoord,geometry);  
  }

  template <class F, class T, class CT=T>
  static int intersectSegmentFacet(const F &facet, 
				   const typename F::Coord pCoord[NDIM], 
				   int dim, typename F::Coord otherCoord,
				   T &intersectionCoord)
  {
    typedef typename F::Coord Coord;

    Coord vCoord[F::NVERT][NDIM];
    for (int n=0;n<F::NVERT;++n)
      {
	typename F::Vertex *v=facet.getVertex(n);
	for (int d=0;d<NDIM;++d)
	  vCoord[n][d]=v->getCoord(d);
      }
  
    return Base::template getSegmentFacetIntersection<Coord,T,CT>
      (vCoord,pCoord,dim,otherCoord,intersectionCoord);  
  }
  

  
   /*
  template <class S>
  static int testCheckOld(const S *simplex, const typename S::Coord pCoord[NDIM])
  {   
    typedef typename S::Coord Coord;
    double vx[3];
    double vy[3];
    double px;
    double py;
    for (int i=0;i<3;++i) vx[i]=simplex->getVertex(i)->getCoord(0);
    for (int i=0;i<3;++i) vy[i]=simplex->getVertex(i)->getCoord(1);
    px=pCoord[0];py=pCoord[1];
    return Base::pnpolyC(3,vx,vy,px,py);
    
    //return Base::test(vCoord,&vCoord[NDIM+1][0]);
  }
  */

  /*
  template <class S, class G>
  static int testCheck(const S *simplex, const typename S::Coord pCoord[NDIM], 
		       const G* geometry)
  {
    typedef typename S::Coord Coord;

    Coord vCoord[S::NVERT][NDIM];
    for (int n=0;n<S::NVERT;++n)
      for (int d=0;d<NDIM;++d)
	vCoord[n][d]=simplex->getVertex(n)->getCoord(d);
    printf("BEFORE\n");
    for (int n=0;n<S::NVERT;++n)
      printf("(%.18lg %.18lg)\n",vCoord[n][0],vCoord[n][1]);
    geometry->template checkCoordsConsistency<Coord,S::NVERT,NDIM>(vCoord,pCoord);
    printf("AFTER\n");
    for (int n=0;n<S::NVERT;++n)
      printf("(%.18lg %.18lg)\n",vCoord[n][0],vCoord[n][1]);
    return Base::test(vCoord,pCoord);
  
  }
  */
};

/** \}*/
#include "../internal/namespace.footer"
#endif
