#ifndef __MPI_EXCHANGE_STRUCT_SIMPLICES_GLOBAL_IDS__
#define __MPI_EXCHANGE_STRUCT_SIMPLICES_GLOBAL_IDS__

#include "../../tools/MPI/mpiCommunication.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  template <class T, class S, bool WITH_INDEX = false>
  struct Mpi_SimplicesGlobalIdsT
  {
  };

  template <class T, class S>
  struct Mpi_SimplicesGlobalIdsT<T, S, false>
  {
    typedef Mpi_SimplicesGlobalIdsT<T, S, false> MyType;
    typedef void *MpiStruct;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef S Simplex;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NVERT = S::NVERT;
    static const int IS_INDEXED = 0;
    static const GlobalIdentityValue empty;

    GlobalIdentityValue vertices[NVERT];
    GlobalIdentityValue gid;

    Mpi_SimplicesGlobalIdsT()
    {
    }

    Mpi_SimplicesGlobalIdsT(Simplex *s, int myRank)
    {
      set(s, myRank);
    }

    bool isEmpty() const
    {
      return false;
    }

    void set(Simplex *s, int myRank)
    {
      for (int i = 0; i < NVERT; i++)
        vertices[i] = s->getVertex(i)->getGlobalIdentity().get();
      gid = s->getGlobalIdentity(myRank).get();
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

      mpiDataType.push_back<GlobalIdentityValue>(+NVERT, OFFSETOF(MyType, vertices[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, gid));

      mpiDataType.commit<MyType>();
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(+NVERT); // the '+' is there on purpose ;)
      disp.push_back(OFFSETOF(MyType,vertices[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,gid));

      int n=sizeof(MyType)-(disp.back()+sizeof(gid));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",
               n,sizeof(MyType),
               (disp.back()+sizeof(gid)));
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

}

#include "../../internal/namespace.footer"
#endif
