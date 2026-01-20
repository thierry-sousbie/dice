#ifndef __MPI_EXCHANGE_STRUCT_SHADOW_SIMPLICES_FROM_GHOST__
#define __MPI_EXCHANGE_STRUCT_SHADOW_SIMPLICES_FROM_GHOST__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{

  // QUERY

  template <class T, class S, bool WITH_INDEX = true>
  struct Mpi_ShadowSimplicesFromGhostQueryT
  {
  };

  template <class T, class S>
  struct Mpi_ShadowSimplicesFromGhostQueryT<T, S, true>
  {
    typedef Mpi_ShadowSimplicesFromGhostQueryT<T, S, true> MyType;
    typedef void *MpiStruct;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef S Simplex;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NNEI = Simplex::NNEI;
    static const GlobalIdentityValue empty; // = GlobalIdentity::empty.get();
    static const int IS_INDEXED = 1;

    unsigned int __mpi_base_cell_index;
    int index;

    Mpi_ShadowSimplicesFromGhostQueryT() : index(-1)
    {
    }

    Mpi_ShadowSimplicesFromGhostQueryT(Simplex *s, unsigned int idx) : index(-1)
    {
      set(s, idx);
    }

    bool isEmpty() const
    {
      return (index < 0);
    }

    unsigned long getBaseCellIndex() const
    {
      return __mpi_base_cell_index;
    }

    void set(Simplex *s, unsigned long idx)
    {
      for (int i = 0; i < NNEI; ++i)
      {
        Simplex *nei = s->getNeighbor(i);
        if ((nei != NULL) && (nei->isShadow()))
        {
          index = i;
          __mpi_base_cell_index = idx;
          // s->template print<LOG_STD_ALL>("ref =");
          break;
        }
      }
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, __mpi_base_cell_index));
      mpiDataType.push_back<int>(1, OFFSETOF(MyType, index));
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

      type.push_back( MPI_Type<int>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,index));

      int n=sizeof(MyType)-(disp.back()+sizeof(index));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",
               n,sizeof(MyType),
               (disp.back()+sizeof(index)));
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
  const typename Mpi_ShadowSimplicesFromGhostQueryT<T, S, true>::GlobalIdentityValue Mpi_ShadowSimplicesFromGhostQueryT<T, S, true>::empty =
      Mpi_ShadowSimplicesFromGhostQueryT<T, S, true>::GlobalIdentity::empty.get();

  // REPLY

  template <class T, class S, bool WITH_INDEX = false>
  struct Mpi_ShadowSimplicesFromGhostReplyT
  {
  };

  template <class T, class S>
  struct Mpi_ShadowSimplicesFromGhostReplyT<T, S, false>
  {
    typedef Mpi_ShadowSimplicesFromGhostReplyT<T, S, false> MyType;
    // typedef void* MpiStruct;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef S Simplex;
    typedef typename S::Vertex Vertex;
    typedef typename T::VertexData VertexData;
    typedef typename T::SimplexData SimplexData;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NNEI = Simplex::NNEI;
    static const int NVERT = Simplex::NVERT;
    static const GlobalIdentityValue empty; // = GlobalIdentity::empty.get();
    static const int IS_INDEXED = 1;

    Coord vertexPos[NDIM_W];
    GlobalIdentityValue vertexGid[NVERT];
    GlobalIdentityValue simplexGid;
    GlobalIdentityValue simplexGeneration;
    unsigned long vertexIndex;
    GlobalIdentityValue generation;
    VertexData vData;
    SimplexData sData;

    template <class L>
    void print()
    {
      if (NDIM == 2)
      {
        glb::console->print<L>("GID=(%d,%d) P=(%g %g ), V=(%d,%d) (%d,%d) (%d,%d) (index=%ld)\n",
                               GlobalIdentity(simplexGid).rank(),
                               GlobalIdentity(simplexGid).id(),
                               vertexPos[0], vertexPos[1],
                               GlobalIdentity(vertexGid[0]).rank(),
                               GlobalIdentity(vertexGid[0]).id(),
                               GlobalIdentity(vertexGid[1]).rank(),
                               GlobalIdentity(vertexGid[1]).id(),
                               GlobalIdentity(vertexGid[2]).rank(),
                               GlobalIdentity(vertexGid[2]).id(),
                               vertexIndex);
      }
      else
      {
        glb::console->print<L>("GID=(%d,%d) P=(%g %g %g), V=(%d,%d) (%d,%d) (%d,%d) (%d,%d) (index=%ld)\n",
                               GlobalIdentity(simplexGid).rank(),
                               GlobalIdentity(simplexGid).id(),
                               vertexPos[0], vertexPos[1], vertexPos[2],
                               GlobalIdentity(vertexGid[0]).rank(),
                               GlobalIdentity(vertexGid[0]).id(),
                               GlobalIdentity(vertexGid[1]).rank(),
                               GlobalIdentity(vertexGid[1]).id(),
                               GlobalIdentity(vertexGid[2]).rank(),
                               GlobalIdentity(vertexGid[2]).id(),
                               GlobalIdentity(vertexGid[3]).rank(),
                               GlobalIdentity(vertexGid[3]).id(),
                               vertexIndex);
      }
    }

    Mpi_ShadowSimplicesFromGhostReplyT() : simplexGid(empty)
    {
    }

    Mpi_ShadowSimplicesFromGhostReplyT(Simplex *ref, int neiIndex, int myRank) : simplexGid(empty)
    {
      set(ref, neiIndex, myRank);
    }

    bool isEmpty() const
    {
      return (simplexGid == empty);
    }

    void set(Simplex *ref, int neiIndex, int myRank)
    {
      Simplex *nei = ref->getNeighbor(neiIndex);
      for (int i = 0; i < NVERT; i++) // NVERT==NNEI
      {
        vertexGid[i] = nei->getVertex(i)->getGlobalIdentity().get();
        if (nei->getNeighbor(i) == ref)
          vertexIndex = i;
      }
      // if (vertexIndex>NNEI)
      // 	{
      // 	  glb::console->print<LOG_ERROR>("NEIGHBOR NOT FOUND.\n");
      // 	  exit(0);
      // 	}
      Vertex *v = nei->getVertex(vertexIndex);
      v->getCoords(vertexPos);
      simplexGid = nei->getGlobalIdentity(myRank).get();
      simplexGeneration = nei->getGeneration().get();
      generation = v->getGeneration().get();
      vData = v->getData();
      sData = nei->getData();
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<Coord>(+NDIM_W, OFFSETOF(MyType, vertexPos[0]));
      mpiDataType.push_back<GlobalIdentityValue>(+NVERT, OFFSETOF(MyType, vertexGid[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, simplexGid));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, simplexGeneration));
      mpiDataType.push_back<unsigned long>(1, OFFSETOF(MyType, vertexIndex));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, generation));

      mpiDataType.push_back(VertexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, vData));
      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, sData));

      mpiDataType.commit<MyType>();
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;

      type.push_back( MPI_Type<Coord>::get() );
      blocklen.push_back(+NDIM_W);
      disp.push_back(OFFSETOF(MyType,vertexPos[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(+NVERT);
      disp.push_back(OFFSETOF(MyType,vertexGid[0]));

      type.push_back( MPI_Type<GlobalIdentityValue>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,simplexGid));

      type.push_back( MPI_Type<unsigned long>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,vertexIndex));

      int n=sizeof(MyType)-(disp.back()+sizeof(vertexIndex));
      if (n)
  {
    glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",
               n,sizeof(MyType),
               (disp.back()+sizeof(vertexIndex)));
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
  const typename Mpi_ShadowSimplicesFromGhostReplyT<T, S, false>::GlobalIdentityValue Mpi_ShadowSimplicesFromGhostReplyT<T, S, false>::empty =
      Mpi_ShadowSimplicesFromGhostReplyT<T, S, false>::GlobalIdentity::empty.get();

} // namespace

#include "../../internal/namespace.footer"
#endif
