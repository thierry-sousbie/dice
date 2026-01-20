#ifndef __MPI_EXCHANGE_STRUCT_COARSEN_SHADOWS__
#define __MPI_EXCHANGE_STRUCT_COARSEN_SHADOWS__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  // unfortunately, Mpi_RefineQueryResultT has to be of type POD as it will be
  // used with MPI, so we cannot inherit and have to specialize the indexed
  // and the non indexed type for each and every MPI_ type ...
  template <class T, class S>
  struct Mpi_CoarsenShadowsT
  {
    typedef Mpi_CoarsenShadowsT<T, S> MyType;
    typedef S Simplex;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty;
    static const int IS_INDEXED = 1;

    unsigned int __mpi_base_cell_index;
    char type;
    char rmVertexIndex;
    char partnerIndex;
    char dummy;
    // GlobalIdentityValue vertices[Simplex::NVERT];
    // GlobalIdentityValue neighbors[Simplex::NNEI];
    // GlobalIdentityValue gid;

    Mpi_CoarsenShadowsT() : type(-1)
    {
    }

    Mpi_CoarsenShadowsT(Simplex *p, unsigned int index, char type_)
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
      type = type_;
      __mpi_base_cell_index = index;
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
    }

    // unsigned long getBaseCellIndex()
    // {
    //   PRINT_SRC_INFO(LOG_ERROR);
    //   glb::console->print<LOG_ERROR>("Wrong type of mpiStruct.\n");
    //   glb::console->print<LOG_ERROR>("Use a structure with base cell index (template param WITH_INDEX=true)\n");
    //   exit(-1);

    //   return -1;
    // }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, __mpi_base_cell_index));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, type));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, rmVertexIndex));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, partnerIndex));
      mpiDataType.push_back<char>(1, OFFSETOF(MyType, dummy));
      mpiDataType.commit<MyType>();
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<unsigned int>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,__mpi_base_cell_index));

      type.push_back( MPI_Type<char>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,type));

      type.push_back( MPI_Type<char>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,rmVertexIndex));

      type.push_back( MPI_Type<char>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,partnerIndex));

      type.push_back( MPI_Type<char>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,dummy));

      const int lastElSize = sizeof(dummy);

      int n=sizeof(MyType)-(disp.back()+lastElSize);
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",n,sizeof(MyType),(disp.back()+lastElSize));
    for (int i=0;i<n;i++)
      {
        type.push_back( MPI_BYTE );
        blocklen.push_back(1);
        disp.push_back(disp.back()+1);
      }
  }

      MPI_Type_create_struct(type.size(), &blocklen[0], &disp[0], &type[0], &result);
      MPI_Type_commit(&result);

      return result;
      */
    }
  };

  template <class T, class S>
  const typename Mpi_CoarsenShadowsT<T, S>::GlobalIdentityValue
      Mpi_CoarsenShadowsT<T, S>::empty = Mpi_CoarsenShadowsT<T, S>::GlobalIdentity::empty.get();

} // namespace

#include "../../internal/namespace.footer"
#endif
