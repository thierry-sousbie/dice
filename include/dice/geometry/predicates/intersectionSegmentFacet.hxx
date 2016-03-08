#ifndef SEGMENT_FACET_INTERSECTION_HXX__
#define SEGMENT_FACET_INTERSECTION_HXX__

#include "filterType.hxx"
#include "internal/pointInSimplexBase_2D.hxx"
#include "internal/pointInSimplexBase_3D.hxx"

/**
 * @file 
 * @brief  Definition of a class to test whether a segment and facet of a simplex intersect
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

namespace predicate
{
  
  /** \addtogroup Geometry
   *   \{
   */

  /**
   * \class SegmentFacetIntersectionT
   * \brief A class to robustly test whether an axis aligned segment intersects a facet of 
   * a simplex. Convention for the SoS is that segment coordinate Pi 
   * along dimension i (0<=i<NDIM) is disturbed as follows : 
   *              Pi -> Pi - pow(eps,NDIM-i)
   * 
   * \tparam ND The number of dimensions
   * \tparam filterType The type of predicate to use, see predicate::filterType::Type 
   *  Default value is predicate::filterType::Adaptive
   * \todo Embed in a higher dimensional space 
   */

  template <int ND, int filterType=predicate::filterType::Adaptive>
  class IntersectionSegmentFacetT : protected internal::PointInSimplexBaseT<ND,filterType>
  {
  protected:
    typedef internal::PointInSimplexBaseT<ND,filterType> Base;
  public:
    static const int NDIM = ND;

    template <class F, class G>
    static int test(const F &facet, 
		    const typename F::Coord pCoord[NDIM], 
		    int dim, typename F::Coord otherCoord,
		    const G* geometry,
		    int coordsAreConsistent=-1)
    {
      typedef char T;
      T dummy;
      return test<F,G,T,T,false>(facet,pCoord,otherCoord,geometry,dummy,coordsAreConsistent);
    }

    /**
     * \param facet the facet of a simplex
     * \param pCoord A pointer to the ND coordinates of an extremity of the segment
     * \param dim the dimension the segment is aligned with
     * \param otherCoord The coordinate of the other extremity of the segment along axis \a dim
     * \param geometry A pointer to the geometry of the mesh
     * \param[out] the interpolated coordinate of the intersection point along \a dim
     * \param coordsAreConsistent For periodic boundaries, set to 1 if coordinates are
     * consistent (i.e. on the same side of the periodic box), 0 if they are not and -1 if 
     * you do not know. This is useless for non periodic boundaries.
     * \tparam F The Facet class
     * \tparam G The class describing geometric properties of the bounding box. See e.g.
     *  GeometricPropertiesT
     * \return true if the the segment and facet intersect
     */
    template <class F, class G, class T, class CT=T, bool INTERP=true>
    static int test(const F &facet, 
		    const typename F::Coord pCoord[NDIM], 
		    int dim, typename F::Coord otherCoord,
		    const G* geometry, T &intersectionCoord,
		    int coordsAreConsistent=-1)
    {
      typedef typename F::Coord Coord;

      Coord vCoord[F::NVERT][NDIM];
      for (int n=0;n<F::NVERT;++n)
	{
	  typename F::Vertex *v=facet.getVertex(n);
	  for (int d=0;d<NDIM;++d)
	    vCoord[n][d]=v->getCoord(d);
	}

      return test<Coord,G,T,CT,INTERP>
	(vCoord,pCoord,dim,otherCoord,geometry,intersectionCoord,coordsAreConsistent);

      /*
      if ((coordsAreConsistent>0)|| 
	  ((coordsAreConsistent<0)&& 
	   geometry->template coordsAreConsistent<Coord,F::NVERT,NDIM>(vCoord,pCoord))
	  )
	return Base::template getSegmentFacetIntersection<Coord,T,CT>
	  (vCoord,pCoord,dim,otherCoord,intersectionCoord);
      else
	return Base::template getSegmentFacetIntersection<Coord,T,G,CT>
	  (vCoord,pCoord,dim,otherCoord,intersectionCoord,geometry);  
      */
    }

    template <class F>
    static int test(const F &facet, 
		    const typename F::Coord pCoord[NDIM], 
		    int dim, typename F::Coord otherCoord)
		    
    {
      typedef char T;
      T dummy;
      return test<F,T,T,false>(facet,pCoord,otherCoord,dummy);
    }
    
    /**
     * \param facet the facet of a simplex
     * \param pCoord A pointer to the ND coordinates of an extremity of the segment
     * \param dim the dimension the segment is aligned with
     * \param otherCoord The coordinate of the other extremity of the segment along axis \a dim
     * \param[out] the interpolated coordinate of the intersection point along \a dim
     * \tparam F The Facet class
     * \return true if the the segment and facet intersect
     */
    template <class F, class T, class CT=T, bool INTERP=true>
    static int test(const F &facet, 
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
      
      return test<Coord,T,CT,INTERP>
	(vCoord,pCoord,dim,otherCoord,intersectionCoord);
      /*
      return Base::template getSegmentFacetIntersection<Coord,T,CT>
	(vCoord,pCoord,dim,otherCoord,intersectionCoord);  
      */
    }

    template <class FCT, class G>
    static int test(const FCT (vCoord)[NDIM][NDIM],
		    const FCT pCoord[NDIM], 
		    int dim, FCT otherCoord,
		    const G* geometry,
		    int coordsAreConsistent=-1)		    
    {
      typedef char T;
      T dummy;
      return test<FCT,G,T,T,false>
	(vCoord,pCoord,dim,otherCoord,geometry,dummy,coordsAreConsistent);
    }

    /**
     * \param vCoord (NDIM) pointers the (NDIM) coordinates of each vertex of the facet
     * \param pCoord A pointer to the ND coordinates of an extremity of the segment
     * \param dim the dimension the segment is aligned with
     * \param otherCoord The coordinate of the other extremity of the segment along axis \a dim
     * \param geometry A pointer to the geometry of the mesh
     * \param[out] the interpolated coordinate of the intersection point along \a dim
     * \param coordsAreConsistent For periodic boundaries, set to 1 if coordinates are
     * consistent (i.e. on the same side of the periodic box), 0 if they are not and -1 if 
     * you do not know. This is useless for non periodic boundaries.
     * \tparam FCT the type of the coordinates
     * \tparam G The class describing geometric properties of the bounding box. See e.g.
     *  GeometricPropertiesT
     * \return true if the the segment and facet intersect
     */
    template <class FCT, class G, class T, class CT=T, bool INTERP=true>
    static int test(const FCT (vCoord)[NDIM][NDIM],
		    const FCT pCoord[NDIM], 
		    int dim, FCT otherCoord,
		    const G* geometry, T &intersectionCoord,
		    int coordsAreConsistent=-1)
    {
      typedef FCT Coord;
      
      if ((coordsAreConsistent>0)|| 
	  ((coordsAreConsistent<0)&& 
	   geometry->template coordsAreConsistent<Coord,NDIM,NDIM>(vCoord,pCoord))
	  )
	return Base::template getSegmentFacetIntersection<Coord,T,CT,INTERP>
	  (vCoord,pCoord,dim,otherCoord,intersectionCoord);
      else
	return Base::template getSegmentFacetIntersection<Coord,T,G,CT,INTERP>
	  (vCoord,pCoord,dim,otherCoord,intersectionCoord,geometry);  
    }

    template <class FCT>
    static int test(const FCT (vCoord)[NDIM][NDIM],
		    const FCT pCoord[NDIM], 
		    int dim, FCT otherCoord)
    {
      typedef char T;
      T dummy;
      return test<FCT,T,T,false>(vCoord,pCoord,dim,otherCoord,dummy);
    }
    
    /**
     * \param facet the facet of a simplex
     * \param vCoord (NDIM) pointers the (NDIM) coordinates of each vertex of the facet
     * \param pCoord A pointer to the ND coordinates of an extremity of the segment
     * \param dim the dimension the segment is aligned with
     * \param otherCoord The coordinate of the other extremity of the segment along axis \a dim
     * \param[out] the interpolated coordinate of the intersection point along \a dim
     * \tparam FCT the type of the coordinates
     * \return true if the the segment and facet intersect
     */
    template <class FCT, class T, class CT=T, bool INTERP=true>
    static int test(const FCT (vCoord)[NDIM][NDIM],
		    const FCT pCoord[NDIM], 
		    int dim, FCT otherCoord,
		    T &intersectionCoord)
    {
      typedef FCT Coord;
  
      return Base::template getSegmentFacetIntersection<Coord,T,CT,INTERP>
	(vCoord,pCoord,dim,otherCoord,intersectionCoord);  
    }
  };

  /** \}*/

}

#include "../../internal/namespace.footer"
#endif
