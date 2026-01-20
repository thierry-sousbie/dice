#ifndef __MPI_EXCHANGE_STRUCT_SHARED_TREE_ROOT__
#define __MPI_EXCHANGE_STRUCT_SHARED_TREE_ROOT__

#include "../../tools/MPI/mpiCommunication.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  template <class T, class R, bool WITH_INDEX = false>
  struct Mpi_SharedTreeRootT
  {
  };

  template <class T, class R>
  struct Mpi_SharedTreeRootT<T, R, false>
  {
    typedef Mpi_SharedTreeRootT<T, R, false> MyType;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NNEI = R::NNEI;

    static const GlobalIdentityValue empty;

    GlobalIdentityValue neighbors[NNEI];
    GlobalIdentityValue globalIdentity;
    long weight;

    Mpi_SharedTreeRootT() : weight(0)
    {
      for (int i = 0; i < NNEI; i++)
        neighbors[i] = empty;
    }

    bool isEmpty() const
    {
      return (weight <= 0);
    }

    unsigned long getBaseCellIndex()
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("Wrong type of mpiStruct.\n");
      glb::console->print<LOG_ERROR>("Use a structure with base cell index (template param WITH_INDEX=true, NOT IMPLEMENTED)\n");
      exit(-1);

      return -1;
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<GlobalIdentityValue>(+NNEI, OFFSETOF(MyType, neighbors[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, globalIdentity));
      mpiDataType.push_back<long>(1, OFFSETOF(MyType, weight));

      mpiDataType.commit<MyType>();
      return mpiDataType;

      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(+NNEI); // the '+' is there on purpose ;)
      disp.push_back(OFFSETOF(MyType,neighbors[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,globalIdentity));

      type.push_back( MPI_Type<long>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,weight));

      int n=sizeof(MyType)-(disp.back()+sizeof(weight));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",n,sizeof(MyType),
               (disp.back()+sizeof(weight)));
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

  template <class T, class R>
  const typename Mpi_SharedTreeRootT<T, R, false>::GlobalIdentityValue Mpi_SharedTreeRootT<T, R, false>::empty =
      Mpi_SharedTreeRootT<T, R, false>::GlobalIdentity::empty.get();

  template <class T, class R, bool WITH_INDEX = false>
  struct SharedTreeRootT : public Base<WITH_INDEX>
  {
    typedef SharedTreeRootT<T, R, WITH_INDEX> MyType;
    typedef Mpi_SharedTreeRootT<T, R, WITH_INDEX> MpiStruct;
    typedef R Root;

    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NNEI = R::NNEI;

    Root *root;

    SharedTreeRootT() : root(NULL)
    {
    }

    void set(Root *curRoot, const MpiStruct &str)
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("DUMMY function: SharedTreeRoot::MpiStructs cannot be converted to SharedTreeRoot directly.\n");
      glb::console->print<LOG_ERROR>("This function should never have been called.\n");
      exit(-1);
    }

    void set(Root *curRoot)
    {
      root = curRoot;
    }

    void setEmpty()
    {
      root = NULL;
    }

    bool isEmpty()
    {
      return (root == NULL);
    }

    MpiStruct getMpiStruct()
    {
      MpiStruct result;
      if (!isEmpty())
      {
        root->getNeighborsGlobalIdentity(result.neighbors);
        result.globalIdentity = root->getGlobalIdentity().get();
        result.weight = root->getWeight();
      }
      return result;
    }
  };

}

#include "../../internal/namespace.footer"
#endif
