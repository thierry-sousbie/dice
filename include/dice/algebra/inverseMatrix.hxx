#ifndef __MY_INVERSEMATRIX_HXX__
#define __MY_INVERSEMATRIX_HXX__

#include "./internal/inverseMatrix_implementation.hxx"

/**
 * @file 
 * @brief  A class to invert small matrixes
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Algebra
 *   \{
 */

/**
 * \class InverseMatrixT
 * \brief A class to invert small matrixes
 *
 * \tparam ND The size of the matrix
 * \tparam T  the element type of the matrix
 * \tparam T2 the element type of the output inverse matrix
 */

template <int ND, class T, class T2>
struct InverseMatrixT
{
  static const int NDIM = ND;
  typedef T Type;
  
  /** \brief Invert the matrix 'mat'
   *  \param mat the matrix to invert
   *  \param[out] out the inverted matrix
   *  \return the determinant of 'mat'
   */
  static inline double 
  compute(const T (&mat)[ND][ND], T2 (&out)[ND][ND])
  {
    return internal::InverseMatrixT<ND,T,T2>::compute(mat,out);    
  }
};

/** \}*/
#include "../internal/namespace.footer"
#endif
