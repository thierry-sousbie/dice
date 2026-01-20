#ifndef __MPI_EXCHANGE_STRUCT_COARSEN_GHOSTS__
#define __MPI_EXCHANGE_STRUCT_COARSEN_GHOSTS__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{

  template <class T, class S>
  struct Mpi_CoarsenGhostsT
  {
    typedef Mpi_CoarsenGhostsT<T, S> MyType;
    typedef S Simplex;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;

    typedef typename T::SimplexData SimplexData;

    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty;
    static const int IS_INDEXED = 1;

    unsigned int __mpi_base_cell_index;
    char type;
    char rmVertexIndex;
    char partnerIndex;
    SimplexData sData;

    Mpi_CoarsenGhostsT() : type(-1)
    {
    }

    Mpi_CoarsenGhostsT(Simplex *p, unsigned int index, char type_)
    {
      set(p, index, type_);
    }

    bool isEmpty() const
    {
      return (type < 0);
    }

    unsigned int getBaseCellIndex() const
    {
      return __mpi_base_cell_index;
    }

    void set(Simplex *s, unsigned int index, char type_)
    {
      __mpi_base_cell_index = index;
      type = type_;

      if (type == 0)
      {
        rmVertexIndex = s->cache.c[2];
        partnerIndex = s->cache.c[3];
      }
      else
      {
        rmVertexIndex = s->cache.c[0];
        partnerIndex = s->cache.c[1];
      }

      sData = s->getData();
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, __mpi_base_cell_index));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, type));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, rmVertexIndex));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, partnerIndex));

      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, sData));

      mpiDataType.commit<MyType>();
      return mpiDataType;
    }
  };

  template <class T, class S>
  const typename Mpi_CoarsenGhostsT<T, S>::GlobalIdentityValue
      Mpi_CoarsenGhostsT<T, S>::empty = Mpi_CoarsenGhostsT<T, S>::GlobalIdentity::empty.get();

} // namespace

#include "../../internal/namespace.footer"
#endif
