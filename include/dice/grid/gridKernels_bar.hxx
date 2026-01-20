#ifndef __GRID_KERNELS_BAR_HXX__
#define __GRID_KERNELS_BAR_HXX__

#include "../tools/helpers/helpers.hxx"
#include "gridKernelsOperations.hxx"
#include "gridKernels_helpers.hxx"

/**
 * @file
 * @brief Definition of generic bar kernel class used to create differentiation and
 * interpolation kernels over a grid
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

namespace gridKernel
{
  /**
   * \class BarStencilBaseT
   * This implements functionalities of a bar stencil but it designed to be inherited from so that one
   * can easily specilize some members. Note that this class is NOT virtual in order to allow template
   * specilization, so be carefull with which member you call ;)
   */
  template <int ND, int SZ, typename ST>
  class BarStencilBaseT
  {
  public:
    typedef BarStencilBaseT<ND, SZ, ST> MyType;
    typedef ST Data;

    static const int SIZE = SZ;
    static const int NDIM = ND;
    static const int SIZE_IS_ODD = SZ & 1;
    static const int SIZE_IS_EVEN = 1 - SIZE_IS_ODD;
    static const int HALF_SIZE_LEFT = (SZ - 1) / 2;
    static const int HALF_SIZE_RIGHT = HALF_SIZE_LEFT + SIZE_IS_EVEN;
    static const int HALF_SIZE_RIGHT_PLUS_ONE = HALF_SIZE_RIGHT + 1;
    static const int NVALUES = SIZE;

  protected:
    Data val[NVALUES];
    long stride[ND];
    int dim[ND];
    size_t lowerLeftShift[NDIM];

    void set(Data *stencil)
    {
      std::copy_n(stencil, NVALUES, val);
    }

    template <class WT>
    void set(Data value, WT *w)
    {
      val[w[0]] = value;
    }

    template <class WT>
    void set(Data value, int i, int j = 0, int k = 0, int l = 0)
    {
      val[i] = value;
    }

  public:
    BarStencilBaseT()
    {
    }

    // NOTE THAT C_ORDER WILL BE TERRIBLE FOR THE CACHE ...
    template <typename IT1, typename IT2, bool FORTRAN_ORDER = true>
    void initialize(IT1 *p_dim, IT2 *p_stride, int strideFactor = 1)
    {
      std::copy_n(p_dim, ND, dim);
      std::copy_n(p_stride, ND, stride);
      for (int i = 0; i < NDIM; ++i)
        stride[i] *= strideFactor;
      if (!FORTRAN_ORDER)
        std::reverse(stride, stride + ND);

      for (int i = 0; i < ND; ++i)
        lowerLeftShift[i] = stride[i] * HALF_SIZE_LEFT;
    }

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t getCorePixelCoords(const DT coords[NDIM], const DT2 origin[NDIM],
                              const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      return gridKernel::getCorePixelCoords<NDIM, PERIODIC, CellCenteredValues, SIZE_IS_ODD>(coords, origin, dim, boxSizeInv, stride, output_pos);
    }

  protected:
    // NOTE: dx is the distance to the central value if SIZE is ODD and the distance to the closest
    // left value if it is EVEN.
    // => if SIZE is ODD  => -0.5 <= dx < 0.5
    // => if SIZE is EVEN =>  0   <= dx < 1
    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t getCorePixelCoordsAndCoefs(const DT coords[NDIM], const DT2 origin[NDIM],
                                      const DT3 boxSizeInv[NDIM], IT output_pos[NDIM],
                                      double dx[NDIM])
    {
      return gridKernel::getCorePixelCoordsAndCoefs<NDIM, PERIODIC, CellCenteredValues, SIZE_IS_ODD>(coords, origin, dim, boxSizeInv, stride, output_pos, dx);
    }

  private:
    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(DT *__restrict ptr, DT2 *value, int dimIndex,
                                 int disp, int len, int nSub,
                                 hlp::ConstantValue<IMPRINT>) const
    {
      DT val = static_cast<DT>(*value);
      for (int i = 0; i < len; ++i, ptr += stride[dimIndex])
        (*ptr) = val;
    }

    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(const DT *__restrict ptr, DT2 *result, int dimIndex,
                                 int disp, int len, int nSub,
                                 hlp::ConstantValue<APPLY>) const
    {
      const Data *__restrict curVal = val + disp;
      if (nSub == 0)
        (*result) = 0;

      for (int i = 0; i < len; ++i, ptr += stride[dimIndex], ++curVal)
        (*result) += (*curVal) * (*ptr);
    }

    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(const DT *__restrict ptr, DT2 *result, int dimIndex,
                                 int disp, int len, int nSub,
                                 hlp::ConstantValue<COPY>) const
    {
      DT2 *__restrict curVal = result + disp;

      for (int i = 0; i < len; ++i, ptr += stride[dimIndex], ++curVal)
        (*curVal) = (*ptr);
    }

  protected:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class DT2, class IT>
    int base_apply_operation_boundary(DT *ptr, DT2 *value, int dimIndex, IT *pos) const
    {
      typedef hlp::ConstantValue<OPERATION> Operation;
      int mask = 0;

      if (pos[dimIndex] < HALF_SIZE_LEFT)
      {
        mask |= (1 << dimIndex);
        double len = HALF_SIZE_LEFT - pos[dimIndex];
        apply_generic_operation<PERIODIC>(ptr - pos[dimIndex] * stride[dimIndex],
                                          value, dimIndex,
                                          len, SIZE - len,
                                          0, Operation());
        if (PERIODIC)
          apply_generic_operation<PERIODIC>(ptr + (dim[dimIndex] - HALF_SIZE_LEFT) * stride[dimIndex],
                                            value, dimIndex,
                                            0, len,
                                            1, Operation());
      }

      if (pos[dimIndex] > (dim[dimIndex] - HALF_SIZE_RIGHT_PLUS_ONE))
      {
        mask |= (1 << (dimIndex + ND));
        double len = dim[dimIndex] - pos[dimIndex] + HALF_SIZE_LEFT;
        apply_generic_operation<PERIODIC>(ptr - pos[dimIndex] * stride[dimIndex],
                                          value, dimIndex,
                                          0, len,
                                          0, Operation());
        if (PERIODIC)
          apply_generic_operation<PERIODIC>(ptr + (dim[dimIndex] - HALF_SIZE_LEFT) * stride[dimIndex],
                                            value, dimIndex,
                                            len, SIZE - len,
                                            1, Operation());
      }

      return mask;
    }

    template <class DT, class DT2>
    void base_apply_operation(DT *__restrict ptr, DT2 *value, int dimIndex,
                              hlp::ConstantValue<IMPRINT>) const
    {
      DT val = static_cast<DT>(*value);
      for (int i = 0; i < SIZE; ++i, ptr += stride[dimIndex])
        (*ptr) = val;
    }

    template <class DT, class DT2>
    void base_apply_operation(const DT *__restrict ptr, DT2 *result, int dimIndex,
                              hlp::ConstantValue<APPLY>) const
    {
      const Data *__restrict curVal = val;

      (*result) = 0;
      for (int i = 0; i < SIZE; ++i, ptr += stride[dimIndex], ++curVal)
        (*result) += (*curVal) * (*ptr);
    }

    template <class DT, class DT2>
    void base_apply_operation(const DT *__restrict ptr, DT2 *result, int dimIndex,
                              hlp::ConstantValue<COPY>) const
    {
      DT2 *__restrict curVal = result;

      for (int i = 0; i < SIZE; ++i, ptr += stride[dimIndex], ++curVal)
        (*curVal) = (*ptr);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void base_apply_operation(DT *data, VT *value, int dimIndex) const
    {
      base_apply_operation(data - lowerLeftShift[dimIndex], value, dimIndex, hlp::ConstantValue<OPERATION>());
    }
  };

  /**
   * \class BarStencilT
   * Generic implementation for any order and dimension. This only uses BarStencilBaseT
   * member functions, and may be specialized for improved performence.
   */
  template <int ND, int SZ, typename ST = double>
  class BarStencilT : public BarStencilBaseT<ND, SZ, ST>
  {
  protected:
    typedef BarStencilBaseT<ND, SZ, ST> B;

  public:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT value, int dimIndex, PT *pos) const
    {
      apply<PERIODIC, OPERATION>(data, &value, dimIndex, pos);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT value, int dimIndex) const
    {
      apply<OPERATION>(data, &value, dimIndex);
    }

    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT *value, int dimIndex, PT *pos) const
    {
      if (!B::base_apply_operation_boundary<PERIODIC, OPERATION>(data, value, dimIndex, pos))
        apply<OPERATION>(data, value, dimIndex);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT *value, int dimIndex) const
    {
      B::template base_apply_operation<OPERATION>(data, value, dimIndex);
    }
  };

} // gridKernel namespace

/** \}*/
#include "../internal/namespace.footer"
#endif
