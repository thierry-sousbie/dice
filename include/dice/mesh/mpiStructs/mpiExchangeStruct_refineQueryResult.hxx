#ifndef __MPI_EXCHANGE_STRUCT_REFINE_QUERY_RESULT__
#define __MPI_EXCHANGE_STRUCT_REFINE_QUERY_RESULT__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  // unfortunately, Mpi_RefineQueryResultT has to be of type POD as it will be
  // used with MPI, so we cannot inherit and have to specialize the indexed
  // and the non indexed type for each and every MPI_ type ...
  template <class T, bool WITH_INDEX = false>
  struct Mpi_RefineQueryResultT
  {
    typedef Mpi_RefineQueryResultT<T, WITH_INDEX> MyType;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;

    GlobalIdentityValue v1;
    GlobalIdentityValue v2;
    double value;

    Mpi_RefineQueryResultT() : value(0)
    {
    }

    bool isEmpty() const
    {
      return (value == 0);
    }

    unsigned long getBaseCellIndex()
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("Wrong type of mpiStruct.\n");
      glb::console->print<LOG_ERROR>("Use a structure with base cell index (template param WITH_INDEX=true)\n");
      exit(-1);

      return -1;
    }

    static MPI_Datatype createMpiStructType()
    {
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back(MPI_Type<GlobalIdentityValue>::get());
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType, v1));

      type.push_back(MPI_Type<GlobalIdentityValue>::get());
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType, v2));

      type.push_back(MPI_Type<double>::get());
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType, value));

      int n = sizeof(MyType) - (disp.back() + sizeof(value));
      if (n)
      {
        glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n", n, sizeof(MyType), (disp.back() + sizeof(value)));
        for (int i = 0; i < n; i++)
        {
          type.push_back(MPI_BYTE);
          blocklen.push_back(1);
          disp.push_back(disp.back() + 1);
        }
      }

      MPI_Type_create_struct(type.size(), &blocklen[0], &disp[0], &type[0], &result);
      MPI_Type_commit(&result);

      return result;
    }
  };

  template <class T>
  struct Mpi_RefineQueryResultT<T, true>
  {
    typedef Mpi_RefineQueryResultT<T, true> MyType;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;

    GlobalIdentityValue v1;
    GlobalIdentityValue v2;
    unsigned int __mpi_base_cell_index;
    float value;

    Mpi_RefineQueryResultT() : value(0)
    {
    }

    bool isEmpty() const
    {
      return (value == 0);
    }

    unsigned long getBaseCellIndex()
    {
      return __mpi_base_cell_index;
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, v1));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, v2));
      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, __mpi_base_cell_index));
      mpiDataType.push_back<float>(1, OFFSETOF(MyType, value));

      mpiDataType.commit<MyType>();
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,v1));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,v2));

      type.push_back( MPI_Type<unsigned int>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,__mpi_base_cell_index));

      type.push_back( MPI_Type<float>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,value));

      int n=sizeof(MyType)-(disp.back()+sizeof(value));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",n,sizeof(MyType),(disp.back()+sizeof(value)));
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

  template <class T, class SH, class S, bool WITH_INDEX = false>
  struct RefineQueryResultT : public Base<WITH_INDEX>
  {
    typedef RefineQueryResultT<T, SH, S, WITH_INDEX> MyType;
    typedef Mpi_RefineQueryResultT<T, WITH_INDEX> MpiStruct;
    typedef SH SegmentHandle;
    typedef typename SegmentHandle::HandledType Segment;
    typedef typename Segment::Vertex Vertex;
    typedef S Simplex;

    double value;
    SegmentHandle handle;

    RefineQueryResultT() : value(0)
    {
    }

    void set(Simplex *curSimplex, const MpiStruct &str)
    {
      value = str.value;
      if (isEmpty())
        return;
      Vertex *v1 = curSimplex->getVertexByGlobalIdentity(str.v1);
      Vertex *v2 = curSimplex->getVertexByGlobalIdentity(str.v2);
      handle = SegmentHandle(Segment(v1, v2, curSimplex));
    }

    void set(double value_, SegmentHandle handle_)
    {
      handle = handle_;
      value = value_;
    }

    void setEmpty()
    {
      value = 0;
    }

    bool isEmpty()
    {
      return (value == 0);
    }

    MpiStruct getMpiStruct()
    {
      MpiStruct result;
      if (!isEmpty())
      {
        Base<WITH_INDEX>::setMpiBaseCellIndex(result);
        result.value = value;
        result.v1 = handle->getVertex(0)->getGlobalIdentity().get();
        result.v2 = handle->getVertex(1)->getGlobalIdentity().get();
      }
      return result;
    }
  };

} // namespace

#include "../../internal/namespace.footer"
#endif
