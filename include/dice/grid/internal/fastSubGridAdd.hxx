#ifndef __FAST_SUBGRID_ADD_HXX__
#define __FAST_SUBGRID_ADD_HXX__

#include "../../internal/namespace.header"

namespace internal
{

  template <int SZ, int NDIM, typename DT>
  class FastFixedSubGridAddT
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM]);
  };

  template <int NDIM, typename DT>
  class FastFixedSubGridAddT<1, NDIM, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      (*data) += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<2, 1, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<2, 2, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<2, 3, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data += stride[2];
      data[0] += value;
      data[1] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<4, 1, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<4, 2, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
    }
  };

  template <typename DT>
  class FastFixedSubGridAddT<4, 3, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[2];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[2];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[2];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
      data += stride[1];
      data[0] += value;
      data[1] += value;
      data[2] += value;
      data[3] += value;
    }
  };

  template <int SZ, typename DT>
  class FastFixedSubGridAddT<SZ, 1, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      for (int i = 0; i < SZ; i += 4)
      {
        data[i + 0] += value;
        data[i + 1] += value;
        data[i + 2] += value;
        data[i + 3] += value;
      }
    }
  };

  template <int SZ, typename DT>
  class FastFixedSubGridAddT<SZ, 2, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      long delta = -stride[1];
      for (int j = 0; j < SZ; ++j)
      {
        delta += stride[1];
        for (int i = delta; i < delta + SZ; i += 4)
        {
          data[i + 0] += value;
          data[i + 1] += value;
          data[i + 2] += value;
          data[i + 3] += value;
        }
      }
    }
  };

  template <int SZ, typename DT>
  class FastFixedSubGridAddT<SZ, 3, DT>
  {
    template <typename IT>
    void add(DT *data, DT value, IT stride[NDIM])
    {
      long delta = -stride[1] - stride[2];
      for (int k = 0; k < SZ; ++k)
      {
        delta += stride[2];
        for (int j = 0; j < SZ; ++j)
        {
          delta += stride[1];
          for (int i = delta; i < delta + SZ; i += 4)
          {
            data[i + 0] += value;
            data[i + 1] += value;
            data[i + 2] += value;
            data[i + 3] += value;
          }
        }
      }
    }
  };

  template <int NDIM, typename DT>
  class FastSubGridAddT
  {
    template <typename IT>
    void add(int level, DT *data, DT value, IT stride[NDIM])
    {
      // if (level==0) return FastFixedSubGridAddT<1,NDIM,DT>::add(data,value,stride);
      switch (level)
      {
      case 0:
        return FastFixedSubGridAddT<1, NDIM, DT>::add(data, value, stride);
        break;
      case 1:
        return FastFixedSubGridAddT<2, NDIM, DT>::add(data, value, stride);
        break;
      case 2:
        return FastFixedSubGridAddT<4, NDIM, DT>::add(data, value, stride);
        break;
      case 3:
        return FastFixedSubGridAddT<8, NDIM, DT>::add(data, value, stride);
        break;
      case 4:
        return FastFixedSubGridAddT<16, NDIM, DT>::add(data, value, stride);
        break;
      case 5:
        return FastFixedSubGridAddT<32, NDIM, DT>::add(data, value, stride);
        break;
      case 6:
        return FastFixedSubGridAddT<64, NDIM, DT>::add(data, value, stride);
        break;
      case 7:
        return FastFixedSubGridAddT<128, NDIM, DT>::add(data, value, stride);
        break;
      case 8:
        return FastFixedSubGridAddT<256, NDIM, DT>::add(data, value, stride);
        break;
      case 9:
        return FastFixedSubGridAddT<512, NDIM, DT>::add(data, value, stride);
        break;
      case 10:
        return FastFixedSubGridAddT<1024, NDIM, DT>::add(data, value, stride);
        break;
      case 11:
        return FastFixedSubGridAddT<2048, NDIM, DT>::add(data, value, stride);
        break;
      case 12:
        return FastFixedSubGridAddT<4096, NDIM, DT>::add(data, value, stride);
        break;
      case 13:
        return FastFixedSubGridAddT<8192, NDIM, DT>::add(data, value, stride);
        break;
      case 14:
        return FastFixedSubGridAddT<16384, NDIM, DT>::add(data, value, stride);
        break;
      case 15:
        return FastFixedSubGridAddT<32768, NDIM, DT>::add(data, value, stride);
        break;
      case 16:
        return FastFixedSubGridAddT<65536, NDIM, DT>::add(data, value, stride);
        break;
      case 17:
        return FastFixedSubGridAddT<131072, NDIM, DT>::add(data, value, stride);
        break;
      case 18:
        return FastFixedSubGridAddT<262144, NDIM, DT>::add(data, value, stride);
        break;
      case default:
        printf("ERROR in internal::FastSubGridAddT : case not handled (N=%ld)\n", 1 << level);
        exit(-1);
      }
    }
  };

} // internal
#include "../../internal/namespace.footer"
#endif
