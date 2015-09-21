#ifndef __ORIENTATION_PREDICATE_HXX__
#define __ORIENTATION_PREDICATE_HXX__

/**
 * @file 
 * @brief Predicates for the position of a point with respect to a codim 1 plane.
 * @author Thierry Sousbie
 */

#include "filterType.hxx"
#include "internal/orientationBase_2D.hxx"
#include "internal/orientationBase_3D.hxx"

#include "../../internal/namespace.header"

namespace predicate
{

  /** \addtogroup Geometry
   *   \{
   */
  
  /**
   * \class OrientationT
   * \brief A class that defines a predicate on the position of a point with respect to 
   * a codim 1 plane.
   * => segment in 2D, plane in 3D, ...
   * Note that the predicate never returns that a point is exactly on the plane. 
   * Instead simulation of simplicity is used to consistently decide which side a point 
   * exactly on the plane should be pushed to. This is achieved by symbolically 
   * perturbing the queried point coordinates Pi as :
   *          Pi -> Pi - pow(eps,NDIM-i)
   *
   * \tparam ND The number of dimension of the embedding space
   * \tparam filterType The type of filter for the predicate 
   * (see predicate::filterType::Type). 
   * Default value is predicate::filterType::Adaptive (i.e. uses exact arithmetic 
   * when necessary only)
   */
  template <int ND,int filterType = predicate::filterType::Adaptive>
  class OrientationT {};
    
  /** \brief 2D specialisation of OrientationT */
  template <int filterType>
  class OrientationT<2,filterType> 
  {
  public:
    static const int NDIM=2;
    
    /** \brief test whether point r is above or below segments [pq] (only for NDIM=2)
     *  \return 0 if below, 1 if above.
     */  
    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM])
    {
      return internal::OrientationBaseT<2,filterType>::
	template test<T,CT>(p,q,r);
    }
	
    /** \brief test whether point r is above or below segments [pq] (only for NDIM=2)
     *  \return 0 if below, 1 if above.
     */  
    template <class T, class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM],const G* geometry)
    {
      return internal::OrientationBaseT<2,filterType>::
	template test<T,G,CT>(p,q,r,geometry);
    }	
  };

  /** \brief 3D specialisation of OrientationT */
  template <int filterType>
  class OrientationT<3,filterType>
  {
  public:
    static const int NDIM=3;

    /** \brief test whether point s is above or below triangle [pqr]
     *  \return 0 if below, 1 if above.
     */  
    template <class T, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM], const T s[NDIM])
    {
      return internal::OrientationBaseT<3,filterType>::
	template test<T,CT>(p,q,r,s);
    }

    /** \brief test whether point s is above or below triangle [pqr]
     *  \return 0 if below, 1 if above.
     */  
    template <class T,class G, class CT=T>
    static int test(const T p[NDIM],const T q[NDIM],const T r[NDIM], const T s[NDIM],
		    const G* geometry)
    {
      return internal::OrientationBaseT<3,filterType>::
	template test<T,G,CT>(p,q,r,s,geometry);
    }
  };

  /** \}*/
 

}


#include "../../internal/namespace.footer"
#endif
