#ifndef __SIMPLEX_INTERPOLATOR_HXX__
#define __SIMPLEX_INTERPOLATOR_HXX__

#include "internal/simplexInterpolatorBase.hxx"

#include "../geometry/barycentricCoordinates.hxx"
#include "../geometry/geometricProperties.hxx"

#include "../tools/helpers/helpers.hxx"

#include "../internal/namespace.header"

template < int ND , class GP = DummyGeometricPropertiesT<ND> >
class SimplexInterpolatorT
{
  typedef internal::SimplexInterpolatorBaseT<ND> Interpolator;

public:
  typedef BarycentricCoordinatesT<ND,GP> BarycentricCoordinates;
  typedef double Coord;
  
  static const int NDIM = ND;
 
  template <class T>
  SimplexInterpolatorT(const T * const coords[NDIM+1],const GP &geometry):
    bc(coords,geometry)
  {}

  template <class T>
  SimplexInterpolatorT(const T coords[NDIM+1][NDIM],const GP &geometry):
    bc(coords,geometry)
  {}
  /*
  template <class S>
  SimplexInterpolatorT(const S* simplex,const GP &geometry):
    bc(simplex,geometry)
  {}
  */
  template <class T>
  SimplexInterpolatorT(const T * const coords[NDIM+1]):
    bc(coords)
  {}

  template <class T>
  SimplexInterpolatorT(const T coords[NDIM+1][NDIM]):
    bc(coords)
  {}
  /*
  template <class S>
  SimplexInterpolatorT(const S* simplex):
    bc(simplex)
  {}
  */
  template <class T, class T2, int order=1>
  T2 interpolate(const T coords[ND], const T2 value[ND+1]) const
  {
    T lambda[NDIM+1];
    bc.computeLambda(coords,lambda);
    return Interpolator::template interpolate<T,T2>(lambda,value,hlp::ConstantValue<order>());
  }

  template <class T, class T2, int order=1>
  T2 interpolateInside(const T coords[ND], const T2 value[ND+1]) const
  {
    T lambda[NDIM+1];
    bc.computeLambda(coords,lambda);
    return Interpolator::template interpolateInside<T,T2,order>(lambda,value);
  }

  /*
  template <class T, class T2>
  T2 interpolate(const T *coords, const T2 *value, int order) const
  {
    T lambda[NDIM+1];
    bc.computeLambda(coords,lambda);
    return Interpolator::interpolate(lambda,value,order);
  }
  */
private:
  const BarycentricCoordinates bc;
};

#include "../internal/namespace.footer"
#endif
