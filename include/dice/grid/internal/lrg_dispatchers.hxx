#ifndef __LOCAL_REGULAR_GRID_DISPATCHERS_HXX__
#define __LOCAL_REGULAR_GRID_DISPATCHERS_HXX__

#include "../../internal/namespace.header"

namespace internal
{
  template <int order, int diff_order = 0>
  struct KernelDispacherT;

  template <>
  struct KernelDispacherT<0, 0>
  {
    template <class M, class CT, bool CellCenteredValues = true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {
      me.applyKernel(me.gridKernel_NGP, coords, result, fieldIndex);
    }
  };

  template <>
  struct KernelDispacherT<1, 0>
  {
    template <class M, class CT, bool CellCenteredValues = true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {
      me.applyKernel(me.gridKernel_CIC, coords, result, fieldIndex);
    }
  };

  template <>
  struct KernelDispacherT<2, 0>
  {
    template <class M, class CT, bool CellCenteredValues = true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {
      me.applyKernel(me.gridKernel_TSC, coords, result, fieldIndex);
    }
  };
  /*
  template <> struct KernelDispacherT<0,0> {
    template <class M, class CT, bool CellCenteredValues=true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {me.applyKernel(me.gridKernel_grad_NGP,coords,result,fieldIndex);}
  };
  */
  template <>
  struct KernelDispacherT<1, 1>
  {
    template <class M, class CT, bool CellCenteredValues = true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {
      me.applyKernel(me.gridKernel_grad_CIC, coords, result, fieldIndex);
    }
  };

  template <>
  struct KernelDispacherT<2, 1>
  {
    template <class M, class CT, bool CellCenteredValues = true>
    static void interpolate(const M &me, CT *coords, double *result, int fieldIndex)
    {
      me.applyKernel(me.gridKernel_grad_TSC, coords, result, fieldIndex);
    }
  };
} // internal

#include "../../internal/namespace.footer"
#endif
