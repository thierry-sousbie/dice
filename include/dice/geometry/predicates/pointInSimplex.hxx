#ifndef __POINT_IN_SIMPLEX_HXX__
#define __POINT_IN_SIMPLEX_HXX__

#include "filterType.hxx"
#include "internal/pointInSimplexBase_2D.hxx"
#include "internal/pointInSimplexBase_3D.hxx"

/**
 * @file 
 * @brief  Definition of a class to test whether a point is inside our outside a simplex
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

namespace predicate
{
  
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

      return test(vCoord,pCoord,geometry,coordsAreConsistent);
      /*
      if ((coordsAreConsistent>0) ||	
	  ((coordsAreConsistent<0)&&
	   geometry->template coordsAreConsistent<Coord,S::NVERT,NDIM>(vCoord,pCoord))
	  )
	return Base::test(vCoord,pCoord);
      else
	return Base::test(vCoord,pCoord,geometry);
      */
    }
  

    /**
     * \brief Test if a point is inside or outside a simplex. Note that the point cannot be 
     *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
     * \param vCoord (NDIM+1) pointers the (NDIM) coordinates of each vertex of the simplex
     * \param pCoord A pointer to the ND coordinates of the point to test
     * \param geometry A pointer to the geometry of the mesh
     * \param coordsAreConsistent For periodic boundaries, set to 1 if coordinates are
     * consistent (i.e. on the same side of the periodic box), 0 if they are not and -1 if 
     * you do not know. This is useless for non periodic boundaries.
     * \tparam T The type of the coordinates
     * \tparam G The class describing geometric properties of the bounding box. See e.g.
     *  GeometricPropertiesT
     * \return true if the point with coords \a pCoord is INSIDE the simplex, false otherwise
     */
    template <class T, class G>
    static int test(const T (vCoord)[NDIM+1][NDIM], const T pCoord[NDIM], 
		    const G* geometry, int coordsAreConsistent=-1)
    {
      if ((coordsAreConsistent>0) ||	
	  ((coordsAreConsistent<0)&&
	   geometry->template coordsAreConsistent<T,NDIM+1,NDIM>(vCoord,pCoord))
	  )
	return Base::test(vCoord,pCoord);
      else
	return Base::test(vCoord,pCoord,geometry);
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
    static int test(const T (vCoord)[NDIM+1][NDIM], const T pCoord[NDIM])
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
    static int test(const T (vCoord)[NDIM+2][NDIM])
    {
      return Base::test(vCoord,&vCoord[NDIM+1][0]);
    }

  };

  /** \}*/

}

#include "../../internal/namespace.footer"
#endif
