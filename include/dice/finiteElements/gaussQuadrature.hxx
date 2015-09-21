#ifndef __GAUSS_QUADRATURE_HXX__
#define __GAUSS_QUADRATURE_HXX__

#include "internal/specializedGQ.hxx"

/**
 * @file 
 * @brief Defines a class for computing gauss quadrature over simplices
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup FiniteElements
 *   \{
 */

/**
 * \class GaussQuadratureT
 * A class to compute gauss quadrature over a simplex given a functor F defining an 
 * operator() that takes NDIM double argument corresponding to the barycentric coordinates
 * and return the value of the function at the given coordinates.
 * \tparam ND  Number of dimensions
 */
template <int ND>
class GaussQuadratureT
{
private:
  //typedef internal::gaussQuadrature::EmptyFunctor EmptyFunctor;

public:
  
  static const int NDIM=ND;
 
  /** \brief Computes the quadrature of the function defined by functor over the simplex.
   *  \param functor A functor whose operator() takes NDIM double valued arguments 
   *  corresponding to the barycentric coordinates and returns the value of the function 
   *  to integrate at that point.
   *  \tparam F The type of the functor, defining a simplexFETypeE::Type static 
   *  constant SIMPLEX_TYPE indicating the type of finite element simplex to 
   *  integrate over. 
   * \tparam DEGREE The degree of the function to integrate. Quadrature will be exact 
   *  if the function is a polynomial of degree <= DEGREE, approximate otherwise.
   * \return the quadrature of the function over the simplex
   * \see simplexFETypeE
   */
  template <class F, int DEGREE=F::QUADRATURE_DEGREE>
  static double compute(F &functor)
  {
    typedef internal::gaussQuadrature::SpecializedGQT<ND,DEGREE,F::FE_TYPE> GQ;
    return GQ::get(functor);
  }

  /** \brief Computes the quadrature over the simplex of the functions defined by 
   *  the functors in the static list functorList .
   *  \param functorList A static list of functors whose operator() takes NDIM double 
   *  valued arguments corresponding to the barycentric coordinates and returns the 
   *  value of the function to integrate at that point.
   *  \param[out] result the quadratures computed for each functor in the list
   *  \tparam FL The type of the static list of functors, each functor defining a 
   *  static constant SIMPLEX_TYPE of type simplexFETypeE::Type indicating the type of 
   *  finite element simplex to integrate over and a static integer value 
   *  QUADRATURE_DEGREE defining the degree of the quadrature. 
   *  The static list can be built using dice::hlp::makeObjectList.
   * \see simplexFETypeE
   */
  template <class FL>
  static void computeList(FL &functorList, double result[FL::SIZE])
  {
    typedef typename hlp::ConstantValue< (FL::SIZE>0) > Continue;
    computeListHelper(functorList,result,Continue());
  }

private:

  template <class FL>
  static void computeListHelper(FL &functorList, double *result, hlp::ConstantValue<false>)
  {}
  
  template <class FL>
  static void computeListHelper(FL &functorList, double *result, hlp::ConstantValue<true>)
  {   
    typedef typename hlp::ConstantValue< (FL::INDEX>0) > Continue;

    (*result)=compute(functorList.getObject());    
    computeListHelper(functorList.getNext(),result+1,Continue());
  }
};

/** \}*/
#include "../internal/namespace.footer"
#endif

