#ifndef __MY_INVERSEMATRIX_IMPL_HXX__
#define __MY_INVERSEMATRIX_IMPL_HXX__

#include "../determinant.hxx"
#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

namespace internal
{

  template <int ND, class T, class T2>
  struct InverseMatrixT;

  template <class T, class T2>
  struct InverseMatrixT<1, T, T2>
  {
    static const int NDIM = 1;
    typedef T Type;
    typedef typename hlp::MostPreciseFP<T, T2>::Result CT;

    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM])
    {
      return compute(mat, out, typename hlp::SameType<T, CT>::Result());
    }

  private:
    static inline CT
    compute(const CT (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsTrue)
    {
      CT det = DeterminantT<NDIM, CT>::compute(mat);
      if (det != 0)
      {
        out[0][0] = CT(1.0) / det;
      }
      return det;
    }
    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsFalse)
    {
      CT mat2[NDIM][NDIM];
      for (int i = 0; i < NDIM; ++i)
        for (int j = 0; j < NDIM; ++j)
          mat2[i][j] = mat[i][j];
      return compute(mat2, out, hlp::IsTrue());
    }
  };

  template <class T, class T2>
  struct InverseMatrixT<2, T, T2>
  {
    static const int NDIM = 2;
    typedef T Type;
    typedef typename hlp::MostPreciseFP<T, T2>::Result CT;

    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM])
    {
      return compute(mat, out, typename hlp::SameType<T, CT>::Result());
    }

  private:
    static inline CT
    compute(const CT (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsTrue)
    {
      CT det = DeterminantT<NDIM, CT>::compute(mat);
      if (det != 0.0)
      {
        CT det_inv = CT(1.0) / det;

        out[0][0] = mat[1][1] * det_inv;
        out[0][1] = -mat[0][1] * det_inv;
        out[1][0] = -mat[1][0] * det_inv;
        out[1][1] = mat[0][0] * det_inv;
      }
      return det;
    }
    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsFalse)
    {
      CT mat2[NDIM][NDIM];
      for (int i = 0; i < NDIM; ++i)
        for (int j = 0; j < NDIM; ++j)
          mat2[i][j] = mat[i][j];
      return compute(mat2, out, hlp::IsTrue());
    }
  };

  template <class T, class T2>
  struct InverseMatrixT<3, T, T2>
  {
    static const int NDIM = 3;
    typedef T Type;
    typedef typename hlp::MostPreciseFP<T, T2>::Result CT;

    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM])
    {
      return compute(mat, out, typename hlp::SameType<T, CT>::Result());
    }

  private:
    static inline CT
    compute(const CT (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsTrue)
    {
      CT det = DeterminantT<NDIM, CT>::compute(mat);
      if (det != 0)
      {
        CT det_inv = CT(1.0) / det;

        out[0][0] = (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) * det_inv;
        out[0][1] = (mat[2][1] * mat[0][2] - mat[0][1] * mat[2][2]) * det_inv;
        out[0][2] = (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]) * det_inv;

        out[1][0] = (mat[2][0] * mat[1][2] - mat[1][0] * mat[2][2]) * det_inv;
        out[1][1] = (mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2]) * det_inv;
        out[1][2] = (mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]) * det_inv;

        out[2][0] = (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]) * det_inv;
        out[2][1] = (mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1]) * det_inv;
        out[2][2] = (mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]) * det_inv;
      }
      return det;
    }
    static inline CT
    compute(const T (&mat)[NDIM][NDIM], T2 (&out)[NDIM][NDIM], hlp::IsFalse)
    {
      CT mat2[NDIM][NDIM];
      for (int i = 0; i < NDIM; ++i)
        for (int j = 0; j < NDIM; ++j)
          mat2[i][j] = mat[i][j];
      return compute(mat2, out, hlp::IsTrue());
    }
  };

}

#include "../../internal/namespace.footer"
#endif
