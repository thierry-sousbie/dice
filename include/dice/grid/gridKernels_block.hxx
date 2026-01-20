#ifndef __GRID_KERNELS_BLOCK_HXX__
#define __GRID_KERNELS_BLOCK_HXX__

#include "../tools/helpers/helpers.hxx"
#include "gridKernelsOperations.hxx"
#include "gridKernels_helpers.hxx"

/**
 * @file
 * @brief Definition of generic block kernel class used to create differentiation and
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
   * \class BlockStencilBaseT
   * This implements all functionalities of a block stencil but it designed to be inherited from so that one
   * can easily specilize some members. Note that this class is NOT virtual in order to allow template
   * specilization, so be carefull with which member you call ;)
   */
  template <int ND, int SZ, typename ST>
  class BlockStencilBaseT
  {
  public:
    typedef BlockStencilBaseT<ND, SZ, ST> MyType;
    typedef ST Data;

    static const int SIZE = SZ;
    static const int NDIM = ND;
    static const int SIZE_IS_ODD = SZ & 1;
    static const int SIZE_IS_EVEN = 1 - SIZE_IS_ODD;
    static const int HALF_SIZE_LEFT = (SZ - 1) / 2;
    static const int HALF_SIZE_RIGHT = HALF_SIZE_LEFT + SIZE_IS_EVEN;
    static const int HALF_SIZE_RIGHT_PLUS_ONE = HALF_SIZE_RIGHT + 1;
    static const int NVALUES = hlp::IntPower<SIZE, NDIM>::Result::value;

    // => divide by 0 for SZ=1
    // static const long LOWER_LEFT_SHIFT = HALF_SIZE_LEFT*(1-hlp::IntPower<SIZE,ND>::Result::value)/(1-SIZE);

  protected:
    Data val[NVALUES];
    long stride[ND];
    int dim[ND];
    size_t lowerLeftShift;

    void set(Data *stencil)
    {
      std::copy_n(stencil, NVALUES, val);
    }

    template <class WT>
    void set(Data value, WT *w)
    {
      long index = 0;
      long sizeStride = 1;
      // if (FORTRAN_ORDER)
      for (int i = 0; i < ND; ++i)
      {
        index += sizeStride * w[i];
        sizeStride *= SIZE;
      }
      // else
      // for (int i=0;i<ND;++i) index += sizeStride[i]*w[ND-i-1];

      val[index] = value;
    }

    template <class WT>
    void set(Data value, int i, int j = 0, int k = 0, int l = 0)
    {
      int w[4] = {i, j, k, l};
      set<WT>(value, w);
    }

  public:
    BlockStencilBaseT()
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

      lowerLeftShift = 0;

      for (int i = 0; i < ND; ++i)
        lowerLeftShift += stride[i] * HALF_SIZE_LEFT;
      /*
      if (FORTRAN_ORDER)
  {
    for (int i=0;i<ND;++i)
      lowerLeftShift+=stride[i]*HALF_SIZE_LEFT;
  }
      else
  {
    for (int i=0;i<ND;++i)
      lowerLeftShift+=stride[ND-i-1]*HALF_SIZE_LEFT;
  }
      */
    }

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t getCorePixelCoords(const DT coords[NDIM], const DT2 origin[NDIM], const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {

      return gridKernel::getCorePixelCoords<NDIM, PERIODIC, CellCenteredValues, SIZE_IS_ODD>(coords, origin, dim, boxSizeInv, stride, output_pos);

      /*
      size_t adress = 0;
      for (int i=0;i<NDIM;++i)
  {
    double tmp = (coords[i]-origin[i])*boxSizeInv[i]*dim[i];
    if (SIZE_IS_EVEN && CellCenteredValues) tmp-=0.5;

    double k;
    if (PERIODIC)
      {
        if (tmp<0) tmp+=dim[i];
        else if (tmp>=dim[i]) tmp-=dim[i];

        modf(tmp,&k);
      }
    else
      {
        modf(tmp,&k);
        if (tmp<0) k=0;
      }
    output_pos[i]=static_cast<IT>(k);
    adress += stride[i]*output_pos[i];
  }
      return adress;
      */
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

      /*
      size_t adress = 0;
      for (int i=0;i<NDIM;++i)
  {
    double tmp = (coords[i]-origin[i])*boxSizeInv[i]*dim[i];
    if (SIZE_IS_EVEN && CellCenteredValues) tmp-=0.5;

    double k;
    if (PERIODIC)
      {
        if (tmp<0) {tmp+=dim[i];}
        else if (tmp>=dim[i]) tmp-=dim[i];

        dx[i]=modf(tmp,&k);
        // dx[i][1] = modf(tmp,&k);
        // dx[i][0] = 1.0 - dx[i][1];
      }
    else
      {
        dx[i]=modf(tmp,&k);
        // dx[i][1] = modf(tmp,&k);
        // dx[i][0] = 1.0 - dx[i][1];
        if (tmp<0)
    {
      k=0;
      dx[i]=0.0;
      // dx[i][0]=1.0;
      // dx[i][1]=0.0;
    }
        if (tmp+1.0 >= dim[i])
    {
      dx[i]=1.0;
      // dx[i][0]=1.0;
      // dx[i][1]=0.0;
    }
      }
    output_pos[i]=static_cast<IT>(k);
    adress += stride[i]*output_pos[i];

    if (SIZE_IS_ODD && CellCenteredValues) dx[i]-=0.5;
  }
      return adress;
      */
    }

  protected:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class DT2, class IT>
    void subDivide(int curDim, DT *__restrict ptr, DT2 *value, IT *pos, int mask,
                   int disp[NDIM], int len[NDIM], int &nSub) const
    {
      typedef hlp::ConstantValue<OPERATION> Operation;

      if (curDim == ND)
        return apply_generic_operation<PERIODIC>(ptr, value, disp, len, nSub++, Operation());

      if (mask & (1 << curDim))
      {
        len[curDim] = HALF_SIZE_RIGHT_PLUS_ONE + pos[curDim];
        disp[curDim] = SIZE - len[curDim];
        subDivide<PERIODIC, OPERATION>(curDim + 1, ptr - pos[curDim] * stride[curDim],
                                       value, pos, mask, disp, len, nSub);

        if (PERIODIC)
        {
          disp[curDim] = 0;
          len[curDim] = SIZE - len[curDim];
          subDivide<PERIODIC, OPERATION>(curDim + 1, ptr + (dim[curDim] - HALF_SIZE_LEFT) * stride[curDim],
                                         value, pos, mask, disp, len, nSub);
        }
      }
      else if (mask & (1 << (curDim + ND)))
      {
        disp[curDim] = 0;
        len[curDim] = dim[curDim] - pos[curDim] + HALF_SIZE_LEFT;
        subDivide<PERIODIC, OPERATION>(curDim + 1, ptr - HALF_SIZE_LEFT * stride[curDim],
                                       value, pos, mask, disp, len, nSub);

        if (PERIODIC)
        {
          disp[curDim] = len[curDim];
          len[curDim] = SIZE - len[curDim];
          subDivide<PERIODIC, OPERATION>(curDim + 1, ptr - pos[curDim] * stride[curDim],
                                         value, pos, mask, disp, len, nSub);
        }
      }
      else
      {
        len[curDim] = SIZE;
        disp[curDim] = 0;
        subDivide<PERIODIC, OPERATION>(curDim + 1, ptr - HALF_SIZE_LEFT * stride[curDim],
                                       value, pos, mask, disp, len, nSub);
      }
    }

    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(DT *__restrict ptr, DT2 *value,
                                 int disp[NDIM], int len[NDIM], int nSub,
                                 hlp::ConstantValue<IMPRINT>) const
    {
      long delta[NDIM - 1];
      DT curVal = static_cast<DT>(*value);

      for (int i = 0; i < ND - 1; ++i)
        delta[i] = stride[i + 1] - (len[i] * stride[i]);

      if (NDIM == 4)
      {
        for (int l = 0; l < len[3]; ++l, ptr += delta[2])
          for (int k = 0; k < len[2]; ++k, ptr += delta[1])
            for (int j = 0; j < len[1]; ++j, ptr += delta[0])
              for (int i = 0; i < len[0]; ++i, ptr += stride[0])
                (*ptr) = curVal;
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < len[2]; ++k, ptr += delta[1])
          for (int j = 0; j < len[1]; ++j, ptr += delta[0])
            for (int i = 0; i < len[0]; ++i, ptr += stride[0])
              (*ptr) = curVal;
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < len[1]; ++j, ptr += delta[0])
          for (int i = 0; i < len[0]; ++i, ptr += stride[0])
            (*ptr) = curVal;
      }
      else
      {
        for (int i = 0; i < len[0]; ++i, ptr += stride[0])
          (*ptr) = curVal;
      }
    }

    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(const DT *__restrict ptr, DT2 *result,
                                 int disp[NDIM], int len[NDIM], int nSub,
                                 hlp::ConstantValue<APPLY>) const
    {
      const Data *__restrict curVal = val;
      long delta[NDIM - 1];
      long deltaKernel[NDIM - 1];
      long kernelStride = 1;

      if (nSub == 0)
        (*result) = 0;

      for (int i = 0; i < ND - 1; ++i)
      {
        delta[i] = stride[i + 1] - (len[i] * stride[i]);
        deltaKernel[i] = kernelStride * (SIZE - len[i]);
        curVal += disp[i] * kernelStride;
        kernelStride *= SIZE;
      }
      curVal += disp[ND - 1] * kernelStride;

      if (NDIM == 4)
      {
        for (int l = 0; l < len[3]; ++l, ptr += delta[2], curVal += deltaKernel[2])
          for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
            for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
              for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
                (*result) += (*curVal) * (*ptr);
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
          for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
            for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
              (*result) += (*curVal) * (*ptr);
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
          for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
            (*result) += (*curVal) * (*ptr);
      }
      else
      {
        for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
          (*result) += (*curVal) * (*ptr);
      }

      // printf("APPLYING (%d %d) / (%d %d)\n %e %e %e %e %e %e %e %e %e\n => %e\n",
      // 	     disp[0],disp[1],len[0],len[1],
      // 	     val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7],val[8],(*result));
    }

    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(const DT *__restrict ptr, DT2 *result,
                                 int disp[NDIM], int len[NDIM], int nSub,
                                 hlp::ConstantValue<COPY>) const
    {
      DT2 *__restrict curVal = result;
      long delta[NDIM - 1];
      long deltaKernel[NDIM - 1];
      long kernelStride = 1;

      for (int i = 0; i < ND - 1; ++i)
      {
        delta[i] = stride[i + 1] - (len[i] * stride[i]);
        deltaKernel[i] = kernelStride * (SIZE - len[i]);
        curVal += disp[i] * kernelStride;
        kernelStride *= SIZE;
      }
      curVal += disp[ND - 1] * kernelStride;

      if (NDIM == 4)
      {
        for (int l = 0; l < len[3]; ++l, ptr += delta[2], curVal += deltaKernel[2])
          for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
            for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
              for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
                (*curVal) = (*ptr);
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
          for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
            for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
              (*curVal) = (*ptr);
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
          for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
            (*curVal) = (*ptr);
      }
      else
      {
        for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
          (*curVal) = (*ptr);
      }

      // printf("APPLYING (%d %d) / (%d %d)\n %e %e %e %e %e %e %e %e %e\n => %e\n",
      // 	     disp[0],disp[1],len[0],len[1],
      // 	     val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7],val[8],(*result));
    }

    // Imprint if non null !!!!
    template <bool PERIODIC, class DT, class DT2>
    void apply_generic_operation(DT *__restrict ptr, DT2 *value,
                                 int disp[NDIM], int len[NDIM], int nSub,
                                 hlp::ConstantValue<INTERNAL1>) const
    {
      const Data *__restrict curVal = val;
      const DT castedVal = static_cast<DT>(*value);
      long delta[NDIM - 1];
      long deltaKernel[NDIM - 1];
      long kernelStride = 1;

      for (int i = 0; i < ND - 1; ++i)
      {
        delta[i] = stride[i + 1] - (len[i] * stride[i]);
        deltaKernel[i] = kernelStride * (SIZE - len[i]);
        curVal += disp[i] * kernelStride;
        kernelStride *= SIZE;
      }
      curVal += disp[ND - 1] * kernelStride;

      if (NDIM == 4)
      {
        for (int l = 0; l < len[3]; ++l, ptr += delta[2], curVal += deltaKernel[2])
          for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
            for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
              for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
                if ((*curVal) != 0)
                  (*ptr) = castedVal;
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < len[2]; ++k, ptr += delta[1], curVal += deltaKernel[1])
          for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
            for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
              if ((*curVal) != 0)
                (*ptr) = castedVal;
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < len[1]; ++j, ptr += delta[0], curVal += deltaKernel[0])
          for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
            if ((*curVal) != 0)
              (*ptr) = castedVal;
      }
      else
      {
        for (int i = 0; i < len[0]; ++i, ptr += stride[0], ++curVal)
          if ((*curVal) != 0)
            (*ptr) = castedVal;
      }

      // printf("APPLYING (%d %d) / (%d %d)\n %e %e %e %e %e %e %e %e %e\n => %e\n",
      // 	     disp[0],disp[1],len[0],len[1],
      // 	     val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7],val[8],(*result));
    }

  protected:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class DT2, class IT>
    int base_apply_operation_boundary(DT *data, DT2 *value, IT *pos) const
    {
      int mask = 0;
      for (int i = 0; i < ND; ++i)
      {
        if (pos[i] < HALF_SIZE_LEFT)
          mask |= (1 << i);
        if (pos[i] > (dim[i] - HALF_SIZE_RIGHT_PLUS_ONE))
          mask |= (1 << (i + ND));
      }

      if (mask)
      {
        int len[NDIM];
        int disp[NDIM];
        int nSub = 0;
        subDivide<PERIODIC, OPERATION>(0, data, value, pos, mask, disp, len, nSub);
      }

      return mask;
    }

    template <KernelOp OPERATION, class DT, class VT>
    void base_apply_operation(DT *data, VT *value) const
    {
      base_apply_operation(data - lowerLeftShift, value, hlp::ConstantValue<OPERATION>());
    }

    template <class DT, class DT2>
    void base_apply_operation(DT *__restrict ptr, DT2 *value, hlp::ConstantValue<IMPRINT>) const
    {
      long delta[NDIM - 1];
      DT curVal = static_cast<DT>(*value);

      for (int i = 0; i < ND - 1; ++i)
        delta[i] = stride[i + 1] - (SIZE * stride[i]);

      if (NDIM == 4)
      {
        for (int l = 0; l < SIZE; ++l, ptr += delta[2])
          for (int k = 0; k < SIZE; ++k, ptr += delta[1])
            for (int j = 0; j < SIZE; ++j, ptr += delta[0])
              for (int i = 0; i < SIZE; ++i, ptr += stride[0])
                (*ptr) = curVal;
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < SIZE; ++k, ptr += delta[1])
          for (int j = 0; j < SIZE; ++j, ptr += delta[0])
            for (int i = 0; i < SIZE; ++i, ptr += stride[0])
              (*ptr) = curVal;
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < SIZE; ++j, ptr += delta[0])
          for (int i = 0; i < SIZE; ++i, ptr += stride[0])
            (*ptr) = curVal;
      }
      else
      {
        for (int i = 0; i < SIZE; ++i, ptr += stride[0])
          (*ptr) = curVal;
      }
    }

    template <class DT, class DT2>
    void base_apply_operation(const DT *__restrict ptr, DT2 *result, hlp::ConstantValue<APPLY>) const
    {
      (*result) = 0;
      const Data *__restrict curVal = val;

      long delta[NDIM - 1];
      for (int i = 0; i < ND - 1; ++i)
        delta[i] = stride[i + 1] - (SIZE * stride[i]);

      if (NDIM == 4)
      {
        for (int l = 0; l < SIZE; ++l, ptr += delta[2])
          for (int k = 0; k < SIZE; ++k, ptr += delta[1])
            for (int j = 0; j < SIZE; ++j, ptr += delta[0])
              for (int i = 0; i < SIZE; ++i, ptr += stride[0])
              {
                (*result) += (*curVal) * (*ptr);
                ++curVal;
              }
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < SIZE; ++k, ptr += delta[1])
          for (int j = 0; j < SIZE; ++j, ptr += delta[0])
            for (int i = 0; i < SIZE; ++i, ptr += stride[0])
            {
              (*result) += (*curVal) * (*ptr);
              ++curVal;
            }
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < SIZE; ++j, ptr += delta[0])
          for (int i = 0; i < SIZE; ++i, ptr += stride[0])
          {
            (*result) += (*curVal) * (*ptr);
            ++curVal;
          }
      }
      else
      {
        for (int i = 0; i < SIZE; ++i, ptr += stride[0])
        {
          (*result) += (*curVal) * (*ptr);
          ++curVal;
        }
      }
    }

    template <class DT, class DT2>
    void base_apply_operation(const DT *__restrict ptr, DT2 *result, hlp::ConstantValue<COPY>) const
    {
      DT2 *__restrict curVal = result;

      long delta[NDIM - 1];
      for (int i = 0; i < ND - 1; ++i)
        delta[i] = stride[i + 1] - (SIZE * stride[i]);

      if (NDIM == 4)
      {
        for (int l = 0; l < SIZE; ++l, ptr += delta[2])
          for (int k = 0; k < SIZE; ++k, ptr += delta[1])
            for (int j = 0; j < SIZE; ++j, ptr += delta[0])
              for (int i = 0; i < SIZE; ++i, ptr += stride[0])
              {
                (*curVal) = (*ptr);
                ++curVal;
              }
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < SIZE; ++k, ptr += delta[1])
          for (int j = 0; j < SIZE; ++j, ptr += delta[0])
            for (int i = 0; i < SIZE; ++i, ptr += stride[0])
            {
              (*curVal) = (*ptr);
              ++curVal;
            }
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < SIZE; ++j, ptr += delta[0])
          for (int i = 0; i < SIZE; ++i, ptr += stride[0])
          {
            (*curVal) = (*ptr);
            ++curVal;
          }
      }
      else
      {
        for (int i = 0; i < SIZE; ++i, ptr += stride[0])
        {
          (*curVal) = (*ptr);
          ++curVal;
        }
      }
    }

    template <class DT, class DT2>
    void base_apply_operation(DT *__restrict ptr, DT2 *value, hlp::ConstantValue<INTERNAL1>) const
    {
      const Data *__restrict curVal = val;
      const DT castedVal = static_cast<DT>(*value);
      long delta[NDIM - 1];
      for (int i = 0; i < ND - 1; ++i)
        delta[i] = stride[i + 1] - (SIZE * stride[i]);

      if (NDIM == 4)
      {
        for (int l = 0; l < SIZE; ++l, ptr += delta[2])
          for (int k = 0; k < SIZE; ++k, ptr += delta[1])
            for (int j = 0; j < SIZE; ++j, ptr += delta[0])
              for (int i = 0; i < SIZE; ++i, ptr += stride[0], ++curVal)
                if ((*curVal) != 0)
                  (*ptr) = castedVal;
      }
      else if (NDIM == 3)
      {
        for (int k = 0; k < SIZE; ++k, ptr += delta[1])
          for (int j = 0; j < SIZE; ++j, ptr += delta[0])
            for (int i = 0; i < SIZE; ++i, ptr += stride[0], ++curVal)
              if ((*curVal) != 0)
                (*ptr) = castedVal;
      }
      else if (NDIM == 2)
      {
        for (int j = 0; j < SIZE; ++j, ptr += delta[0])
          for (int i = 0; i < SIZE; ++i, ptr += stride[0], ++curVal)
            if ((*curVal) != 0)
              (*ptr) = castedVal;
      }
      else
      {
        for (int i = 0; i < SIZE; ++i, ptr += stride[0], ++curVal)
          if ((*curVal) != 0)
            (*ptr) = castedVal;
      }
    }
  };

  /**
   * \class BlockStencilT
   * Generic implementation for any order and dimension. This only uses BlockStencilBaseT
   * member functions, and should be specialized for improved performences ...
   */
  template <int ND, int SZ, typename ST = double>
  class BlockStencilT : public BlockStencilBaseT<ND, SZ, ST>
  {
  protected:
    typedef BlockStencilBaseT<ND, SZ, ST> B;

  public:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT value, PT *pos) const
    {
      apply<PERIODIC, OPERATION>(data, &value, pos);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT value) const
    {
      apply<OPERATION>(data, &value);
    }

    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT *value, PT *pos) const
    {
      if (B::template base_apply_operation_boundary<PERIODIC, OPERATION>(data, value, pos) == 0)
        apply<OPERATION>(data, value);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT *value) const
    {
      B::template base_apply_operation<OPERATION>(data, value);
    }
  };

  template <typename ST>
  class BlockStencilT<2, 3, ST> : public BlockStencilBaseT<2, 3, ST>
  {
  protected:
    typedef BlockStencilBaseT<2, 3, ST> B;

  public:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT value, PT *pos) const
    {
      apply<PERIODIC, OPERATION>(data, &value, pos);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT value) const
    {
      apply<OPERATION>(data, &value);
    }

    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT *value, PT *pos) const
    {
      if (B::template base_apply_operation_boundary<PERIODIC, OPERATION>(data, value, pos) == 0)
        apply<OPERATION>(data, value);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT *value) const
    {
      applyOperation(data, value, hlp::ConstantValue<OPERATION>());
    }

  protected:
    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<APPLY>) const
    {
      const DT *__restrict d = data - B::lowerLeftShift;
      long s0 = B::stride[0];
      long s0x2 = s0 + s0;

      (*value) = d[0] * B::val[0] + d[s0] * B::val[1] + d[s0x2] * B::val[2];
      d += B::stride[1];
      (*value) += d[0] * B::val[3] + d[s0] * B::val[4] + d[s0x2] * B::val[5];
      d += B::stride[1];
      (*value) += d[0] * B::val[6] + d[s0] * B::val[7] + d[s0x2] * B::val[8];
    }

    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<COPY>) const
    {
      const DT *__restrict d = data - B::lowerLeftShift;
      long s0 = B::stride[0];
      long s0x2 = s0 + s0;

      value[0] = d[0];
      value[1] = d[s0];
      value[2] = d[s0x2];
      d += B::stride[1];
      value[3] = d[0];
      value[4] = d[s0];
      value[5] = d[s0x2];
      d += B::stride[1];
      value[6] = d[0];
      value[7] = d[s0];
      value[8] = d[s0x2];
    }

    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<IMPRINT>) const
    {
      DT *__restrict d = data - B::lowerLeftShift;
      DT v = static_cast<DT>(*value);
      long s0 = B::stride[0];
      long s0x2 = s0 + s0;

      d[0] = v;
      d[s0] = v;
      d[s0x2] = v;
      d += B::stride[1];
      d[0] = v;
      d[s0] = v;
      d[s0x2] = v;
      d += B::stride[1];
      d[0] = v;
      d[s0] = v;
      d[s0x2] = v;
    }
  };

  template <typename ST>
  class BlockStencilT<2, 2, ST> : public BlockStencilBaseT<2, 2, ST>
  {
  protected:
    typedef BlockStencilBaseT<2, 2, ST> B;

  public:
    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT value, PT *pos) const
    {
      apply<PERIODIC, OPERATION>(data, &value, pos);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT value) const
    {
      apply<OPERATION>(data, &value);
    }

    template <bool PERIODIC, KernelOp OPERATION, class DT, class VT, class PT>
    void apply(DT *data, VT *value, PT *pos) const
    {
      if (!B::template base_apply_operation_boundary<PERIODIC, OPERATION>(data, value, pos))
        apply<OPERATION>(data, value);
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT *value) const
    {
      applyOperation(data, value, hlp::ConstantValue<OPERATION>());
    }

  protected:
    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<APPLY>) const
    {
      const DT *__restrict d = data - B::lowerLeftShift;
      long s0 = B::stride[0];

      (*value) = d[0] * B::val[0] + d[s0] * B::val[1];
      d += B::stride[1];
      (*value) += d[0] * B::val[2] + d[s0] * B::val[3];
    }

    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<COPY>) const
    {
      const DT *__restrict d = data - B::lowerLeftShift;
      long s0 = B::stride[0];

      value[0] = d[0];
      value[1] = d[s0];
      d += B::stride[1];
      value[2] = d[0];
      value[3] = d[s0];
    }

    template <class DT, class VT>
    void applyOperation(DT *data, VT *value, hlp::ConstantValue<IMPRINT>) const
    {
      DT *__restrict d = data - B::lowerLeftShift;
      DT v = static_cast<DT>(*value);
      long s0 = B::stride[0];

      d[0] = v;
      d[s0] = v;
      d += B::stride[1];
      d[0] = v;
      d[s0] = v;
    }
  };

} // gridKernel namespace

/** \}*/
#include "../internal/namespace.footer"
#endif
