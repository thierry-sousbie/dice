#ifndef __MPI_EXCHANGE_STRUCT_REPART_SIMPLICES_WITH_CACHE__
#define __MPI_EXCHANGE_STRUCT_REPART_SIMPLICES_WITH_CACHE__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  // unfortunately, Mpi_RefineQueryResultT has to be of type POD as it will be
  // used with MPI, so we cannot inherit and have to specialize the indexed
  // and the non indexed type for each and every MPI_ type ...
  template <class T, class S, bool WITH_INDEX = false>
  struct Mpi_RepartSimplicesWithCacheT
  {
    typedef Mpi_RepartSimplicesWithCacheT<T, S, WITH_INDEX> MyType;
    typedef S Simplex;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef typename T::SimplexData SimplexData;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty;

    GlobalIdentityValue vertices[Simplex::NVERT];
    GlobalIdentityValue neighbors[Simplex::NNEI];
    GlobalIdentityValue gid;
    GlobalIdentityValue generation;
    long cache;
    SimplexData sData;

    Mpi_RepartSimplicesWithCacheT()
    {
    }

    Mpi_RepartSimplicesWithCacheT(Simplex *s, int myRank)
    {
      set(s, myRank);
    }

    bool isEmpty() const
    {
      return false;
    }

    void set(Simplex *s, int myRank)
    {
      // printf("s=%ld\n",(unsigned long)s);
      for (int i = 0; i < Simplex::NVERT; ++i)
        vertices[i] = s->getVertex(i)->getGlobalIdentity().get();
      for (int i = 0; i < Simplex::NNEI; ++i)
      {
        if (s->getNeighbor(i) == NULL)
          neighbors[i] = empty;
        else
          neighbors[i] = s->getNeighbor(i)->getGlobalIdentity(myRank).get();
      }
      gid = s->getGlobalIdentity(myRank).get();
      generation = s->getGeneration().get();
      cache = s->cache.l;
      sData = s->getData();
    }

    unsigned long getBaseCellIndex()
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("Wrong type of mpiStruct.\n");
      glb::console->print<LOG_ERROR>("Use a structure with base cell index (template param WITH_INDEX=true)\n");
      exit(-1);

      return -1;
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;
      // Don't remove the '+' !!!!
      mpiDataType.push_back<GlobalIdentityValue>(+Simplex::NVERT, OFFSETOF(MyType, vertices[0]));
      mpiDataType.push_back<GlobalIdentityValue>(+Simplex::NNEI, OFFSETOF(MyType, neighbors[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, gid));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, generation));
      mpiDataType.push_back<long>(1, OFFSETOF(MyType, cache));

      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, sData));

      mpiDataType.commit<MyType>();
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(+Simplex::NVERT);
      disp.push_back(OFFSETOF(MyType,vertices[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(+Simplex::NNEI);
      disp.push_back(OFFSETOF(MyType,neighbors[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,gid));

      type.push_back( MPI_Type<long>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,cache));

      int n=sizeof(MyType)-(disp.back()+sizeof(cache));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",n,sizeof(MyType),(disp.back()+sizeof(cache)));
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

  template <class T, class S, bool WITH_INDEX>
  const typename Mpi_RepartSimplicesWithCacheT<T, S, WITH_INDEX>::GlobalIdentityValue Mpi_RepartSimplicesWithCacheT<T, S, WITH_INDEX>::empty =
      Mpi_RepartSimplicesWithCacheT<T, S, WITH_INDEX>::GlobalIdentity::empty.get();

} // namespace

#include "../../internal/namespace.footer"
#endif
