#ifndef __BARYCENTRIC_COORDS_TO_SHAPE_COORDS_HXX__
#define __BARYCENTRIC_COORDS_TO_SHAPE_COORDS_HXX__

#include "../../internal/namespace.header"
namespace internal
{
  template <int NDIM, int ORDER>
  struct QuadraticSimplexGetIndex;
  template <int NDIM, int ORDER>
  struct BarycentricCoordsToShapeCoords;

  template <int NDIM>
  struct QuadraticSimplexGetIndex<NDIM, 2>
  {
    static int get(int vi, int si)
    {
      if (si < 0)
        return vi;
      else if (vi > si)
        std::swap(vi, si);

      return NDIM * vi - (vi * (vi - 1)) / 2 + si - vi - 1;
    }
  };

  template <int NDIM>
  struct BarycentricCoordsToShapeCoords<NDIM, 2>
  {
    template <class I, class O>
    static void convert(I in, O out)
    {
      for (int i = 0; i < NDIM + 1; ++i)
        out[i] = in[i] * (in[i] * 2.0 - 1.0);

      int n = NDIM + 1;
      for (int i = 0; i < NDIM; ++i)
        for (int j = i + 1; j < NDIM + 1; ++j, ++n)
          out[n] = in[i] * in[j] * 4.0;
    }
  };

  template <>
  struct BarycentricCoordsToShapeCoords<3, 2>
  {
    template <class I, class O>
    static void convert(I in, O out)
    {
      out[0] = in[0] * (in[0] * 2.0 - 1.0);
      out[1] = in[1] * (in[1] * 2.0 - 1.0);
      out[2] = in[2] * (in[2] * 2.0 - 1.0);
      out[3] = in[3] * (in[3] * 2.0 - 1.0);
      out[4] = in[0] * in[1] * 4.0;
      out[5] = in[0] * in[2] * 4.0;
      out[6] = in[0] * in[3] * 4.0;
      out[7] = in[1] * in[2] * 4.0;
      out[8] = in[1] * in[3] * 4.0;
      out[9] = in[2] * in[3] * 4.0;
    }
  };

  template <>
  struct BarycentricCoordsToShapeCoords<2, 2>
  {
    template <class I, class O>
    static void convert(I in, O out)
    {
      out[0] = in[0] * (in[0] * 2.0 - 1.0);
      out[1] = in[1] * (in[1] * 2.0 - 1.0);
      out[2] = in[2] * (in[2] * 2.0 - 1.0);
      out[3] = in[0] * in[1] * 4.0;
      out[4] = in[0] * in[2] * 4.0;
      out[5] = in[1] * in[2] * 4.0;
    }
  };

}
#include "../../internal/namespace.footer"
#endif
