#ifndef __SIMPLEX_VOLUME_HXX__
#define __SIMPLEX_VOLUME_HXX__

#include "./internal/simplexVolume_implementation.hxx"
#include <cmath>

/**
 * @file
 * @brief  Simplex volume computation
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Geometry
 *   \{
 */

/**
 * \class SimplexVolumeT
 * \brief N-Simplex volume computation
 *
 * \tparam ND  Number of the dimensions of the simplex
 * \tparam NDW Number of dimensions of the embedding space
 * \tparam T   Data type
 */

template <int ND, int NDW, class T>
struct SimplexVolumeT : protected internal::SimplexVolumeT<ND, NDW, T>
{
  typedef internal::SimplexVolumeT<ND, NDW, T> Base;
  typedef SimplexVolumeT<ND, NDW, T> MyType;
  static const int NDIM = Base::NDIM;
  static const int NDIM_W = Base::NDIM_W;

  /** \brief Compute the oriented (i.e. signed) squared volume.
   *  \param base The NDIM base vectors of the simplex (V_i-V_0) where V_i are the
   *  coordinates of the ith Vertex of the oriented simplex.
   *  \return the signed squared volume of the oriented simplex
   */
  static inline T compute_2S(const T (&base)[NDIM][NDIM_W])
  {
    return Base::compute_2S(base);
  }

  /** \brief Compute the non-oriented (i.e. positive) squared volume.
   *  \param base The NDIM base vectors of the simplex (V_i-V_0) where V_i are the
   *  coordinates of the ith Vertex of the simplex.
   *  \return the squared volume of the simplex.
   */
  static inline T compute_2(const T (&base)[NDIM][NDIM_W])
  {
    T tmp = Base::compute_2S(base);
    if (tmp < 0)
      return -tmp;
    else
      return tmp;
    // return std::fabs(Base::compute_2S(base));
  }

  /** \brief Compute the oriented (i.e. signed) volume.
   *  \param base The NDIM base vectors of the simplex (V_i-V_0) where V_i are the
   *  coordinates of the ith Vertex of the oriented simplex.
   *  \return the signed volume of the simplex.
   */
  static inline T compute_S(const T (&base)[NDIM][NDIM_W])
  {
    T result = Base::compute_2S(base);
    return (result < 0) ? (-sqrt(-result)) : sqrt(result);
  }

  /** \brief Compute the non-oriented (i.e. positive) volume.
   *  \param base The NDIM base vectors of the simplex (V_i-V_0) where V_i are the
   *  coordinates of the ith Vertex of the simplex.
   *  \return the volume of the simplex.
   */
  static inline T compute(const T (&base)[NDIM][NDIM_W])
  {
    T tmp = Base::compute_2S(base);
    if (tmp < 0)
      return sqrt(-tmp);
    else
      return sqrt(tmp);
    // return sqrt(std::fabs(Base::compute_2S(base)));
  }
};

template <int ND, class T>
struct SimplexVolumeT<ND, ND, T> : protected internal::SimplexVolumeT<ND, ND, T>
{
  typedef internal::SimplexVolumeT<ND, ND, T> Base;
  typedef SimplexVolumeT<ND, ND, T> MyType;
  static const int NDIM = Base::NDIM;
  static const int NDIM_W = Base::NDIM_W;

  //! Compute the oriented (i.e. signed) squared volume
  static inline T compute_2S(const T (&base)[NDIM][NDIM_W])
  {
    T result = Base::compute_S(base);
    return (result < 0) ? (-result * result) : (result * result);
  }

  //! Compute the non-oriented (i.e. positive) squared volume
  static inline T compute_2(const T (&base)[NDIM][NDIM_W])
  {
    T result = Base::compute_S(base);
    return result * result;
  }

  //! Compute the oriented (i.e. signed) volume
  static inline T compute_S(const T (&base)[NDIM][NDIM_W])
  {
    return Base::compute_S(base);
  }

  //! Compute the non-oriented (i.e. positive) volume
  static inline T compute(const T (&base)[NDIM][NDIM_W])
  {
    T tmp = Base::compute_S(base);
    if (tmp < 0)
      return -tmp;
    else
      return tmp;
    // return std::fabs(Base::compute_S(base));
  }
};

/** \}*/
#include "../internal/namespace.footer"
#endif
