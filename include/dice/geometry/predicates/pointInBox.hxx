#ifndef __POINT_IN_BOX_HXX__
#define __POINT_IN_BOX_HXX__

#include "filterType.hxx"
#include "./internal/pointInBoxBase.hxx"

/**
 * @file 
 * @brief  Definition of a class to test whether a point is inside our outside a box
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

namespace predicate
{
  
  /** \addtogroup Geometry
   *   \{
   */

  /**
   * \class PointInBoxT
   * \brief A class to test whether a point is inside our outside a box
   * 
   * \tparam ND The number of dimensions
   * \tparam filterType The type of predicate to use, see predicate::filterType::Type 
   *  Default value is predicate::filterType::Adaptive
   * \todo Embed in a higher dimensional space 
   */

  template <int ND, int filterType=predicate::filterType::Adaptive>
  class PointInBoxT : protected internal::PointInBoxBaseT<ND,filterType>
  {
  protected:
    typedef internal::PointInBoxBaseT<ND,filterType> Base;
  public:
    static const int NDIM = ND;

    /**
     * \brief Test if a point is inside or outside a box. Note that the point cannot be 
     *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
     * \param box an array of the NDIM coordinates of 2 opposite vertices of the box
     * \param pCoord A pointer to the ND coordinates of the point to test
     * \param geometry A pointer to the geometry of the mesh
     * \param coordsAreConsistent For periodic boundaries, set to 1 if coordinates are
     * consistent (i.e. on the same side of the periodic box), 0 if they are not and -1 if 
     * you do not know. This is useless for non periodic boundaries.
     * \tparam G The class describing geometric properties of the bounding box. See e.g.
     *  GeometricPropertiesT
     * \return true if the point with coords \a pCoord is INSIDE the box, false otherwise
     */
    template <class T, class G>
    static int test(const T box[2][NDIM], const T pCoord[NDIM], 
		    const G* geometry, int coordsAreConsistent=-1)
    {
      if (coordsAreConsistent>0)
	return Base::test(box,pCoord);
      else
	return Base::test(box,pCoord,geometry);
    }
  
 
    /**
     * \brief Test if a point is inside or outside a box. Note that the point cannot be 
     *  on a boundary, in which case a side is arbitrarily (but consistently) chosen.
     * \param box an array of the NDIM coordinates of 2 opposite vertices of the box
     * \param pCoord A pointer to the ND coordinates of the point to test
     * \tparam T The type of the coordinates
     * \return True if the point with coords \a pCoord is INSIDE the box, false otherwise
     */
    template <class T>
    static int test(const T box[2][NDIM], const T pCoord[NDIM])
    {
      return Base::test(box,pCoord);
    }

  };

  /** \}*/

}

#include "../../internal/namespace.footer"
#endif
