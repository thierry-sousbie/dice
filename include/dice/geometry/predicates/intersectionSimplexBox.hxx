#ifndef INTERSECTION_SIMPLEX_BOX_HXX__
#define INTERSECTION_SIMPLEX_BOX_HXX__

#include "filterType.hxx"
#include "./internal/intersectionSimplexBoxBase.hxx"

/**
 * @file 
 * @brief  Definition of a class to test whether a box and a simplex intersect. The class also
 * allows for the retrieval of the intersection points.
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

namespace predicate
{
  
  /** \addtogroup Geometry
   *   \{
   */

  /**
   * \class intersectionSimplexBoxT
   * \brief A class to test whether a box and a simplex intersect. The class also
   * allows for the retrieval of the intersection points.
   * 
   * \tparam ND The number of dimensions
   * \tparam filterType The type of predicate to use, see predicate::filterType::Type 
   *  Default value is predicate::filterType::Adaptive
   * \todo Embed in a higher dimensional space 
   */

  template <int ND, int filterType=predicate::filterType::Adaptive>
  class IntersectionSimplexBoxT : protected internal::IntersectionSimplexBoxBaseT<ND,filterType>
  {
  protected:
    typedef internal::IntersectionSimplexBoxBaseT<ND,filterType> Base;
    
  public:
    static const int NDIM = ND;

    template <class S, class G>
    static int test(const S *simplex, const typename S::Coord box[2][NDIM], 
		    const G* geometry, int coordsAreConsistent=-1)
    {
      typedef typename S::Coord Coord;

      Coord vCoord[S::NVERT][NDIM];
      for (int n=0;n<S::NVERT;++n)
	for (int d=0;d<NDIM;++d)
	  vCoord[n][d]=simplex->getVertex(n)->getCoord(d);

      return test(vCoord,box,geometry,coordsAreConsistent);
    }
    

    template <class T, class G>
    static int test(const T (&vCoord)[NDIM+1][NDIM], const T box[2][NDIM], 
		    const G* geometry, int coordsAreConsistent=-1)
    {
      if (coordsAreConsistent)
	return Base::test(vCoord,box);

      int consistent=0;
      consistent = geometry->template coordsAreConsistent<T,NDIM+1,NDIM>(vCoord,box[0]);
      consistent += geometry->template coordsAreConsistent<T,NDIM+1,NDIM>(vCoord,box[1]);
      
      if (consistent<2) return Base::test(vCoord,box,geometry);
      
      return Base::test(vCoord,box);
    }
  
    template <class T>
    static int test(const T (&vCoord)[NDIM+1][NDIM], const T box[2][NDIM])
    {
      return Base::test(vCoord,box);
    }

  };

  /** \}*/

}

#include "../../internal/namespace.footer"
#endif
