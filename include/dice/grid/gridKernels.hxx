#ifndef __GRID_KERNELS_HXX__
#define __GRID_KERNELS_HXX__

#include "gridKernels_block.hxx"
#include "gridKernels_bar.hxx"

/**
 * @file
 * @brief Definition of interpolation and differentiation kernels over a regular grid
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

/** \namespace gridKernel
 * \brief contains classes for interpolations and finite differences over a regular grid
 */
namespace gridKernel
{

  /**
   * \class NGP
   * \brief Nearest grid point interpolation kernel (order 0)
   */
  template <int ND, typename ST = double>
  class NGP : public BlockStencilT<ND, 1, ST>
  {
  private:
    typedef BlockStencilT<ND, 1, ST> B;

  public:
    static const int NDIM = ND;
    static const int OUT_COUNT = 1;

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t setCoefs(const DT coords[NDIM], const DT2 origin[NDIM],
                    const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      double dx[NDIM];
      size_t address = B::template getCorePixelCoordsAndCoefs<PERIODIC, CellCenteredValues>(coords, origin, boxSizeInv, output_pos, dx);
      B::val[0] = 1;
      return address;
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM])
    {
      B::val[0] = 1;
    }

    template <class DT>
    void setSpacing(DT w[NDIM])
    {
    }

    void setSpacing(double w, int i = -1)
    {
    }
  };

  /**
   * \class CIC
   * \brief Cloud in cell interpolation kernel (order 1)
   */
  template <int ND, typename ST = double>
  class CIC : public BlockStencilT<ND, 2, ST>
  {
  private:
    typedef BlockStencilT<ND, 2, ST> B;

  public:
    static const int NDIM = ND;
    static const int OUT_COUNT = 1;

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t setCoefs(const DT coords[NDIM], const DT2 origin[NDIM],
                    const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      double dx[NDIM];
      size_t address = B::template getCorePixelCoordsAndCoefs<PERIODIC, CellCenteredValues>(coords, origin, boxSizeInv, output_pos, dx);

      setCoefs(dx, hlp::ConstantValue<ND>());
      return address;
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM])
    {
      setCoefs(dx, hlp::ConstantValue<ND>());
    }

    template <class DT>
    void setSpacing(DT w[NDIM])
    {
    }

    void setSpacing(double w, int i = -1)
    {
    }

  private:
    template <class DT>
    void setCoefs(const DT dxR[NDIM], hlp::ConstantValue<1>)
    {
      double dxL = 1.0 - dxR[0];
      B::val[0] = dxR[0];
      B::val[1] = dxL;
    }

    template <class DT>
    void setCoefs(const DT dxR[NDIM], hlp::ConstantValue<2>)
    {
      double dxL[NDIM] = {1.0 - dxR[0], 1.0 - dxR[1]};
      B::val[0] = dxL[0] * dxL[1];
      B::val[1] = dxR[0] * dxL[1];
      B::val[2] = dxL[0] * dxR[1];
      B::val[3] = dxR[0] * dxR[1];
    }

    template <class DT>
    void setCoefs(const DT dxR[NDIM], hlp::ConstantValue<3>)
    {
      double dxL[NDIM] = {1.0 - dxR[0], 1.0 - dxR[1], 1.0 - dxR[2]};

      double du0 = dxL[0] * dxL[1];
      double du1 = dxR[0] * dxL[1];
      double du2 = dxL[0] * dxR[1];
      double du3 = dxR[0] * dxR[1];

      B::val[0] = du0 * dxL[2];
      B::val[1] = du1 * dxL[2];
      B::val[2] = du2 * dxL[2];
      B::val[3] = du3 * dxL[2];
      B::val[4] = du0 * dxR[2];
      B::val[5] = du1 * dxR[2];
      B::val[6] = du2 * dxR[2];
      B::val[7] = du3 * dxR[2];
    }
  };

  /**
   * \class TSC
   * \brief Trianglular shape cloud interpolation kernel (order 2)
   */
  template <int ND, typename ST = double>
  class TSC : public BlockStencilT<ND, 3, ST>
  {
  private:
    typedef BlockStencilT<ND, 3, ST> B;

  public:
    static const int NDIM = ND;
    static const int OUT_COUNT = 1;

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t setCoefs(const DT coords[NDIM], const DT2 origin[NDIM],
                    const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      double dx[NDIM];
      size_t address = B::template getCorePixelCoordsAndCoefs<PERIODIC, CellCenteredValues>(coords, origin, boxSizeInv, output_pos, dx);

      setCoefs(dx, hlp::ConstantValue<ND>());
      // if (debug) printf("dx=[%e %e], %e %e %e %e %e %e %e %e %e\n",dx[0],dx[1],
      // 			B::val[0],B::val[1],B::val[2],
      // 			B::val[3],B::val[4],B::val[5],
      // 			B::val[6],B::val[7],B::val[8]);
      return address;
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM])
    {
      setCoefs(dx, hlp::ConstantValue<ND>());
    }

    template <class DT>
    void setSpacing(DT w[NDIM])
    {
    }

    void setSpacing(double w, int i = -1)
    {
    }

  private:
    template <class DT>
    void setCoefs(const DT dx[NDIM], hlp::ConstantValue<1>)
    {
      double dxL = 0.5 * (0.5 - dx[0]) * (0.5 - dx[0]); // Left point coef
      double dxC = 0.75 - dx[0] * dx[0];                // Central point coef
      double dxR = 0.5 * (0.5 + dx[0]) * (0.5 + dx[0]); // Right point coef

      B::val[0] = dxL;
      B::val[1] = dxC;
      B::val[2] = dxR;
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM], hlp::ConstantValue<2>)
    {
      double dxL[NDIM] = {0.5 * (0.5 - dx[0]) * (0.5 - dx[0]),
                          0.5 * (0.5 - dx[1]) * (0.5 - dx[1])};
      double dxC[NDIM] = {0.75 - dx[0] * dx[0],
                          0.75 - dx[1] * dx[1]};
      double dxR[NDIM] = {0.5 * (0.5 + dx[0]) * (0.5 + dx[0]),
                          0.5 * (0.5 + dx[1]) * (0.5 + dx[1])};

      B::val[0] = dxL[0] * dxL[1];
      B::val[1] = dxC[0] * dxL[1];
      B::val[2] = dxR[0] * dxL[1];

      B::val[3] = dxL[0] * dxC[1];
      B::val[4] = dxC[0] * dxC[1];
      B::val[5] = dxR[0] * dxC[1];

      B::val[6] = dxL[0] * dxR[1];
      B::val[7] = dxC[0] * dxR[1];
      B::val[8] = dxR[0] * dxR[1];
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM], hlp::ConstantValue<3>)
    {
      double dxL[NDIM] = {0.5 * (0.5 - dx[0]) * (0.5 - dx[0]),
                          0.5 * (0.5 - dx[1]) * (0.5 - dx[1]),
                          0.5 * (0.5 - dx[2]) * (0.5 - dx[2])};
      double dxC[NDIM] = {0.75 - dx[0] * dx[0],
                          0.75 - dx[1] * dx[1],
                          0.75 - dx[2] * dx[2]};
      double dxR[NDIM] = {0.5 * (0.5 + dx[0]) * (0.5 + dx[0]),
                          0.5 * (0.5 + dx[1]) * (0.5 + dx[1]),
                          0.5 * (0.5 + dx[2]) * (0.5 + dx[2])};

      double LL = dxL[0] * dxL[1];
      double CL = dxC[0] * dxL[1];
      double RL = dxR[0] * dxL[1];

      double LC = dxL[0] * dxC[1];
      double CC = dxC[0] * dxC[1];
      double RC = dxR[0] * dxC[1];

      double LR = dxL[0] * dxR[1];
      double CR = dxC[0] * dxR[1];
      double RR = dxR[0] * dxR[1];

      B::val[0] = LL * dxL[2];
      B::val[1] = CL * dxL[2];
      B::val[2] = RL * dxL[2];
      B::val[3] = LC * dxL[2];
      B::val[4] = CC * dxL[2];
      B::val[5] = RC * dxL[2];
      B::val[6] = LR * dxL[2];
      B::val[7] = CR * dxL[2];
      B::val[8] = RR * dxL[2];

      B::val[9] = LL * dxC[2];
      B::val[10] = CL * dxC[2];
      B::val[11] = RL * dxC[2];
      B::val[12] = LC * dxC[2];
      B::val[13] = CC * dxC[2];
      B::val[14] = RC * dxC[2];
      B::val[15] = LR * dxC[2];
      B::val[16] = CR * dxC[2];
      B::val[17] = RR * dxC[2];

      B::val[18] = LL * dxR[2];
      B::val[19] = CL * dxR[2];
      B::val[20] = RL * dxR[2];
      B::val[21] = LC * dxR[2];
      B::val[22] = CC * dxR[2];
      B::val[23] = RC * dxR[2];
      B::val[24] = LR * dxR[2];
      B::val[25] = CR * dxR[2];
      B::val[26] = RR * dxR[2];
    }
  };

  /**
   * \class InterpolateT
   * \brief grid interpolator at order ORDER
   * \tparam ND number of dimensions
   * \tparam ORDER order of the interpolator (0=>NGP, 1=>CIC, 2=>TSC, ...) (only implemented up to 2)
   */
  template <int ND, int ORDER, typename ST = double>
  class InterpolateT;

  template <int ND, typename ST>
  class InterpolateT<ND, 0, ST> : public NGP<ND, ST>
  {
  };
  template <int ND, typename ST>
  class InterpolateT<ND, 1, ST> : public CIC<ND, ST>
  {
  };
  template <int ND, typename ST>
  class InterpolateT<ND, 2, ST> : public TSC<ND, ST>
  {
  };

  /**
   * \class CentralDiffT
   * \brief (2*ORDER+1) points central finite differences
   * \tparam ORDER finite difference order : kernel size is (2*ORDER+1)
   * \tparam DIFF_ORDER order of the derivative (only ORDER=1 is implemented)
   */
  template <int ND, int ORDER, int DIFF_ORDER = 1, typename ST = double>
  class CentralDiffT;

  template <int ND, int ORDER, typename ST>
  class CentralDiffT<ND, ORDER, 1, ST> : public BarStencilT<ND, 2 * ORDER + 1, ST>
  {
  private:
    typedef BarStencilT<ND, 2 * ORDER + 1, ST> B;

  public:
    static const int NDIM = ND;
    static const int OUT_COUNT = 1;

    CentralDiffT(double w = 1)
    {
      setCoeffs(hlp::ConstantValue<ORDER>(), w);
    }
    /*
    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t setCoefs(const DT coords[NDIM],const DT2 origin[NDIM],
        const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      double dx[NDIM];
      return B::template getCorePixelCoordsAndCoefs<PERIODIC,CellCenteredValues>
  (coords,origin,boxSizeInv,output_pos,dx);
    }
    */
    void setSpacing(double w)
    {
      setCoeffs(hlp::ConstantValue<ORDER>(), w);
    }

  private:
    void multCoefs(double w_inv)
    {
      for (int i = 0; i < B::SIZE; ++i)
        B::val[i] *= w_inv;
    }

    void setCoeffs(hlp::ConstantValue<0>, double w = 1)
    {
      B::val[0] = 0;
    }

    void setCoeffs(hlp::ConstantValue<1>, double w = 1)
    {
      B::val[0] = -1.0 / 2.0;
      B::val[1] = 0;
      B::val[2] = 1.0 / 2.0;
      if (w != 1)
        multCoefs(1.0 / w);
    }

    void setCoeffs(hlp::ConstantValue<2>, double w = 1)
    {
      B::val[0] = 1.0 / 12.0;
      B::val[1] = -8.0 / 12.0;
      B::val[2] = 0;
      B::val[3] = 8.0 / 12.0;
      B::val[4] = -1.0 / 12.0;
      if (w != 1)
        multCoefs(1.0 / w);
    }

    void setCoeffs(hlp::ConstantValue<3>, double w = 1)
    {
      B::val[0] = -1.0 / 60.0;
      B::val[1] = 9.0 / 60.0;
      B::val[2] = -45.0 / 60.0;
      B::val[3] = 0;
      B::val[4] = 45.0 / 60.0;
      B::val[5] = -9.0 / 60.0;
      B::val[6] = 1.0 / 60.0;
      if (w != 1)
        multCoefs(1.0 / w);
    }

    void setCoeffs(hlp::ConstantValue<4>, double w = 1)
    {
      B::val[0] = 3.0 / 840.0;
      B::val[1] = -32.0 / 840.0;
      B::val[2] = 168.0 / 840.0;
      B::val[3] = -672.0 / 840.0;
      B::val[4] = 0;
      B::val[5] = 672.0 / 840.0;
      B::val[6] = -168.0 / 840.0;
      B::val[7] = 32.0 / 840.0;
      B::val[8] = -3.0 / 840.0;
      if (w != 1)
        multCoefs(1.0 / w);
    }
  };

  template <int ND, int ORDER, int DIFF_ORDER = 1, typename ST = double>
  class InterpolateCentralDiffT : protected BlockStencilT<ND, 3 * ORDER + 1, char>
  {
  private:
    typedef InterpolateT<ND, ORDER, ST> Interpolate;
    // typedef CentralDiffT<1,ORDER,DIFF_ORDER,ST> CentralDiff;
    typedef CentralDiffT<ND, ORDER, DIFF_ORDER, ST> CentralDiff;

    typedef BlockStencilT<ND, 3 * ORDER + 1, char> B;
    typedef ST Data;
    typedef typename B::Data BaseData;

  public:
    static const int NDIM = ND;
    static const int OUT_COUNT = ND;

    InterpolateCentralDiffT()
    {
      long stride = 1;

      // m_center=B::val;
      for (int i = 0; i < NDIM; ++i)
      {
        // m_center += B::HALF_SIZE_LEFT*stride;
        interpDelta[i] = stride * (B::SIZE - Interpolate::SIZE);
        stride *= B::SIZE;
      }

      // int j=0;
      int w[NDIM];
      int dim[NDIM];
      for (int i = 0; i < NDIM; ++i)
      {
        w[i] = 0;
        dim[i] = B::SIZE;
      }
      for (int i = 0; i < B::NVALUES; ++i)
      {
        int count = 0;
        for (int j = 0; j < NDIM; ++j)
          if ((w[j] < ORDER) || (w[j] >= B::SIZE - ORDER))
            count++;

        if (count > 1)
          B::val[i] = 0;
        else
          B::val[i] = 1;

        hlp::getNext<NDIM>(w, dim);
      }
      // std::fill_n(B::val,B::NVALUES,1.0);
    }
    /*
    template <typename IT1, typename IT2, bool FORTRAN_ORDER=true>
    void initialize(IT1 *p_dim, IT2 *p_stride, int strideFactor=1)
    {
      B::initialize(p_dim,p_stride,strideFactor);
      //smallBlock.initialize(p_dim,p_stride,strideFactor);
      long i_dim[NDIM];
      long i_stride[NDIM+1];
      i_stride[0]=1;
      for (int i=0;i<NDIM;++i)
  {
    i_dim[i]=B::SIZE;
    i_stride[i+1]=i_stride[i]*B::SIZE;
  }
      interp.initialize(i_dim,i_stride);

      int c_dim=CentralDiff::SIZE;
      int c_stride=1;
      for (int i=0;i<NDIM;++i)
  centralDiff[i].initialize(&c_dim,&c_stride);
    }
    */
    template <typename IT1, typename IT2, bool FORTRAN_ORDER = true>
    void initialize(IT1 *p_dim, IT2 *p_stride, int strideFactor = 1)
    {
      B::initialize(p_dim, p_stride, strideFactor);

      long i_dim[NDIM];
      long i_stride[NDIM + 1];
      i_stride[0] = 1;
      for (int i = 0; i < NDIM; ++i)
      {
        i_dim[i] = Interpolate::SIZE;
        i_stride[i + 1] = i_stride[i] * Interpolate::SIZE;
      }
      interp.initialize(i_dim, i_stride);

      i_stride[0] = 1;
      for (int i = 0; i < NDIM; ++i)
      {
        i_dim[i] = B::SIZE;
        i_stride[i + 1] = i_stride[i] * B::SIZE;
      }
      for (int i = 0; i < NDIM; ++i)
        centralDiff[i].initialize(i_dim, i_stride);
    }

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t getCorePixelCoords(const DT coords[NDIM], const DT2 origin[NDIM], const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      return B::template getCorePixelCoords<PERIODIC, CellCenteredValues>(coords, origin, boxSizeInv, output_pos);
    }

    template <bool PERIODIC, bool CellCenteredValues, class DT, class DT2, class DT3, class IT>
    size_t setCoefs(const DT coords[NDIM], const DT2 origin[NDIM],
                    const DT3 boxSizeInv[NDIM], IT output_pos[NDIM])
    {
      double dx[NDIM];
      size_t address = B::template getCorePixelCoordsAndCoefs<PERIODIC, CellCenteredValues>(coords, origin, boxSizeInv, output_pos, dx);
      interp.setCoefs(dx);
      return address;
    }

    template <class DT>
    void setCoefs(const DT dx[NDIM])
    {
      interp.setCoefs(dx);
    }

    template <class DT>
    void setSpacing(DT w[NDIM])
    {
      for (int i = 0; i < NDIM; ++i)
        centralDiff[i].setSpacing(w[i]);
    }

    void setSpacing(double w, int dim = -1)
    {
      if (dim < 0)
      {
        for (int i = 0; i < NDIM; ++i)
          centralDiff[i].setSpacing(w);
      }
      else
        centralDiff[dim].setSpacing(w);
    }

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
      my_apply<PERIODIC>(data, value, pos, hlp::ConstantValue<OPERATION>());
    }

    template <KernelOp OPERATION, class DT, class VT>
    void apply(DT *data, VT *value) const
    {
      my_apply(data, value, hlp::ConstantValue<OPERATION>());
    }

  private:
    CentralDiff centralDiff[NDIM];
    Interpolate interp;
    int interpDelta[NDIM];
    // BaseData *m_center;
    /*
    template <class DT, class VT>
    void apply_APPLY(DT *data, VT *value) const
    {
      double interpValue[CentralDiff::SIZE];
      long stride=1;

      for (int i=0;i<NDIM;++i)
  {
    DT * __restrict ptr=data + B::lowerLeftShift - stride*CentralDiff::HALF_SIZE_LEFT;
    for (int j=0;j<CentralDiff::SIZE;++j,ptr+=stride)
      interp.apply<APPLY>(ptr,&interpValue[j]);

    centralDiff[i].apply<APPLY>(interpValue+CentralDiff::HALF_SIZE_LEFT,&value[i],0);
    stride *= B::SIZE;
  }
    }
    */
    template <class DT, class VT>
    void apply_APPLY(DT *data, VT *value) const
    {
      static const long deltaDiff = (B::HALF_SIZE_LEFT - Interpolate::HALF_SIZE_LEFT) *
                                    (1 - hlp::IntPower<B::SIZE, +ND>::Result::value) / (1 - B::SIZE);

      static const long deltaInterp = Interpolate::HALF_SIZE_LEFT *
                                      (1 - hlp::IntPower<Interpolate::SIZE, +ND>::Result::value) / (1 - Interpolate::SIZE);

      double buffer[Interpolate::NVALUES];
      int dim[NDIM];
      std::fill_n(dim, NDIM, +Interpolate::SIZE);

      for (int d = 0; d < NDIM; ++d)
      {
        int w[NDIM] = {0};
        int k = deltaDiff;
        for (int i = 0; i < Interpolate::NVALUES; ++i, ++k)
        {
          centralDiff[d].template apply<APPLY>(&data[k], &buffer[i], d);
          hlp::getNext<NDIM>(w, k, dim, interpDelta);
        }

        interp.template apply<APPLY>(buffer + deltaInterp, &value[d]);
      }
    }

    template <bool PERIODIC, class DT, class VT, class PT>
    void my_apply(DT *data, VT *value, PT *pos, hlp::ConstantValue<APPLY>) const
    {
      Data tmpArr[B::NVALUES];
      B::template apply<PERIODIC, COPY>(data, tmpArr, pos);
      apply_APPLY(tmpArr, value);
    }

    template <class DT, class VT>
    void my_apply(DT *data, VT *value, hlp::ConstantValue<APPLY>) const
    {
      Data tmpArr[B::NVALUES];
      B::template apply<COPY>(data, tmpArr);
      apply_APPLY(tmpArr, value);
    }

    template <bool PERIODIC, class DT, class VT, class PT>
    void my_apply(DT *data, VT *value, PT *pos, hlp::ConstantValue<COPY>) const
    {
      B::template apply<PERIODIC, COPY>(data, value, pos);
    }

    template <class DT, class VT>
    void my_apply(DT *data, VT *value, hlp::ConstantValue<COPY>) const
    {
      B::template apply<COPY>(data, value);
    }

    template <bool PERIODIC, class DT, class VT, class PT>
    void my_apply(DT *data, VT *value, PT *pos, hlp::ConstantValue<IMPRINT>) const
    {
      B::template apply<PERIODIC, INTERNAL1>(data, value, pos);
    }

    template <class DT, class VT>
    void my_apply(DT *data, VT *value, hlp::ConstantValue<IMPRINT>) const
    {
      B::template apply<INTERNAL1>(data, value);
    }
  };

}

/** \}*/
#include "../internal/namespace.footer"
#endif
