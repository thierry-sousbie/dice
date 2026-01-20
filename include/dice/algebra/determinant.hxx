#ifndef __MY_DETERMINANT_HXX__
#define __MY_DETERMINANT_HXX__

#include "./internal/determinant_implementation.hxx"

/**
 * @file
 * @brief  A class to compute the determinant of small matrices
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Algebra
 *   \{
 */

/**
 * \class DeterminantT
 * \brief A class to compute the determinant of small matrices
 *
 * \tparam ND The size of the matrix
 * \tparam T  the element type of the matrix
 */

template <int ND, class T>
struct DeterminantT
{
  static const int NDIM = ND;
  typedef T Type;

  /** \brief Compute the determinant of 'mat'
   *  \param mat the matrix
   *  \return the determinant of 'mat'
   */
  static inline T
  compute(const T (&mat)[ND][ND])
  {
    return internal::DeterminantT<ND, T>::compute(mat);
  }
};

/** \}*/
#include "../internal/namespace.footer"
#endif
