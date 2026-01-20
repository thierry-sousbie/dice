#ifndef __INTERNAL_INTERPOLATE_CIC_HXX__
#define __INTERNAL_INTERPOLATE_CIC_HXX__

#include "../../internal/namespace.header"

namespace internal
{

  template <class DT, int NDIM>
  struct InterpolateCIC;

  template <class DT>
  struct InterpolateCIC<DT, 1>
  {
    static const int NDIM = 1;
    static double compute(const DT *d, const long delta[NDIM], const double dx[NDIM][2])
    {
      double result = d[0] * dx[0][0];
      result += d[delta[0]] * dx[0][1];
      return result;
    }
  };

  template <class DT>
  struct InterpolateCIC<DT, 2>
  {
    static const int NDIM = 2;
    static double compute(const DT *d, const long delta[NDIM], const double dx[NDIM][2])
    {
      double result = d[0] * dx[0][0] * dx[1][0];
      result += d[delta[0]] * dx[0][1] * dx[1][0];
      result += d[delta[1]] * dx[0][0] * dx[1][1];
      result += d[delta[0] + delta[1]] * dx[0][1] * dx[1][1];
      return result;
    }

    template <int Order = 2>
    static void computeGrad(const DT *d, const long delta[NDIM][2],
                            const double dx[NDIM][2], double result[NDIM])
    {
      return computeGrad(d, delta, dx, result, hlp::ConstantValue<Order>());
    }
    // 1/12 -2/3 0 2/3 -1/12
    static void computeGrad(const DT *d, const long delta[NDIM][2],
                            const double dx[NDIM][2], double result[NDIM],
                            hlp::ConstantValue<1>)
    {
      double ddx[4] = {dx[0][0] * dx[1][0],
                       dx[0][1] * dx[1][0],
                       dx[0][0] * dx[1][1],
                       dx[0][1] * dx[1][1]};

      for (int dim = 0; dim < NDIM; ++dim)
      {
        const DT *ptr = d - delta[dim][0];
        result[dim] = (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[0];

        ptr = d + delta[0][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[1];

        ptr = d + delta[1][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[2];

        ptr = d + delta[0][1] + delta[1][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[3];
      }
    }
  };

  template <class DT>
  struct InterpolateCIC<DT, 3>
  {
    static const int NDIM = 3;
    static double compute(const DT *d, const long delta[NDIM], const double dx[NDIM][2])
    {
      double du0 = dx[0][0] * dx[1][0];
      double du1 = dx[0][1] * dx[1][0];
      double du2 = dx[0][0] * dx[1][1];
      double du3 = dx[0][1] * dx[1][1];

      double result = d[0] * du0 * dx[2][0];
      result += d[delta[0]] * du1 * dx[2][0];
      result += d[delta[1]] * du2 * dx[2][0];
      result += d[delta[0] + delta[1]] * du3 * dx[2][0];

      d += delta[2];
      result += d[0] * du0 * dx[2][1];
      result += d[delta[0]] * du1 * dx[2][1];
      result += d[delta[1]] * du2 * dx[2][1];
      result += d[delta[0] + delta[1]] * du3 * dx[2][1];
      return result;
    }

    template <int Order = 2>
    static void computeGrad(const DT *d, const long delta[NDIM][2],
                            const double dx[NDIM][2], double result[NDIM])
    {
      return computeGrad(d, delta, dx, result, hlp::ConstantValue<Order>());
    }

    static void computeGrad(const DT *d, const long delta[NDIM][2],
                            const double dx[NDIM][2], double result[NDIM],
                            hlp::ConstantValue<2>)
    {
      double du0 = dx[0][0] * dx[1][0];
      double du1 = dx[0][1] * dx[1][0];
      double du2 = dx[0][0] * dx[1][1];
      double du3 = dx[0][1] * dx[1][1];

      double ddx[8] = {du0 * dx[2][0], du1 * dx[2][0], du2 * dx[2][0], du3 * dx[2][0],
                       du0 * dx[2][1], du1 * dx[2][1], du2 * dx[2][1], du3 * dx[2][1]};

      for (int dim = 0; dim < NDIM; ++dim)
      {
        const DT *ptr = d - delta[dim][0];
        result[dim] = (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[0];

        ptr = d + delta[0][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[1];

        ptr = d + delta[1][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[2];

        ptr = d + delta[0][1] + delta[1][1] - delta[dim][0];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[3];

        // +delta[2][1]
        ptr = d - delta[dim][0] + delta[2][1];
        result[dim] = (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[4];

        ptr = d + delta[0][1] - delta[dim][0] + delta[2][1];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[5];

        ptr = d + delta[1][1] - delta[dim][0] + delta[2][1];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[6];

        ptr = d + delta[0][1] + delta[1][1] - delta[dim][0] + delta[2][1];
        result[dim] += (ptr[0] - ptr[delta[dim][0]] + ptr[delta[dim][0] + delta[dim][1]]) * ddx[7];
      }
    }
  };

  template <class DT>
  struct InterpolateCIC<DT, 4>
  {
    static const int NDIM = 4;
    static double compute(const DT *d, const long delta[NDIM], const double dx[NDIM][2])
    {
      double du0 = dx[0][0] * dx[1][0];
      double du1 = dx[0][1] * dx[1][0];
      double du2 = dx[0][0] * dx[1][1];
      double du3 = dx[0][1] * dx[1][1];

      double dv0 = dx[2][0] * dx[3][0];
      double dv1 = dx[2][0] * dx[3][1];
      double dv2 = dx[2][1] * dx[3][0];
      double dv3 = dx[2][1] * dx[3][1];

      double result = d[0] * du0 * dv0;
      result += d[delta[0]] * du1 * dv0;
      result += d[delta[1]] * du2 * dv0;
      result += d[delta[0] + delta[1]] * du3 * dv0;

      d += delta[2];
      result += d[0] * du0 * dv1;
      result += d[delta[0]] * du1 * dv1;
      result += d[delta[1]] * du2 * dv1;
      result += d[delta[0] + delta[1]] * du3 * dv1;

      d += delta[3] - delta[2];
      result = d[0] * du0 * dv2;
      result += d[delta[0]] * du1 * dv2;
      result += d[delta[1]] * du2 * dv2;
      result += d[delta[0] + delta[1]] * du3 * dv2;

      d += delta[2];
      result += d[0] * du0 * dv3;
      result += d[delta[0]] * du1 * dv3;
      result += d[delta[1]] * du2 * dv3;
      result += d[delta[0] + delta[1]] * du3 * dv3;
      return result;
    }
  };

} // internal
#include "../../internal/namespace.footer"
#endif
