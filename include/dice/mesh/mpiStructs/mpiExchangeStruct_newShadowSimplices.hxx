#ifndef __MPI_EXCHANGE_STRUCT_NEW_SHADOW_SIMPLICES__
#define __MPI_EXCHANGE_STRUCT_NEW_SHADOW_SIMPLICES__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  template <class T, class S>
  struct Mpi_NewShadowSimplicesT
  {
    typedef Mpi_NewShadowSimplicesT<T, S> MyType;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;

    typedef S Simplex;
    typedef typename S::Vertex Vertex;

    typedef typename T::SimplexData SimplexData;
    typedef typename T::VertexData VertexData;

    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const int NNEI = Simplex::NNEI;
    static const int NVERT = Simplex::NVERT;
    static const GlobalIdentityValue empty;          // = GlobalIdentity::empty.get();
    static const GlobalIdentityValue neighborIsVoid; // = GlobalIdentity::empty.get();
    static const int IS_INDEXED = 1;

    Coord vertexPos[NDIM_W];
    GlobalIdentityValue vertexGid[NVERT];
    GlobalIdentityValue simplexGid;
    GlobalIdentityValue simplexGeneration;
    GlobalIdentityValue generation;
    unsigned int vertexIndex;
    unsigned int __mpi_base_cell_index;
    SimplexData sData;
    VertexData vData;

    unsigned int getBaseCellIndex() const
    {
      return __mpi_base_cell_index;
    }

    template <class L>
    void print()
    {
      glb::console->print<L>("GID=(%d,%d) P=(%g %g %g), V=(%d,%d) (%d,%d) (%d,%d) (index=%d)\n",
                             GlobalIdentity(simplexGid).rank(),
                             GlobalIdentity(simplexGid).id(),
                             vertexPos[0], vertexPos[1], vertexPos[2],
                             GlobalIdentity(vertexGid[0]).rank(),
                             GlobalIdentity(vertexGid[0]).id(),
                             GlobalIdentity(vertexGid[1]).rank(),
                             GlobalIdentity(vertexGid[1]).id(),
                             GlobalIdentity(vertexGid[2]).rank(),
                             GlobalIdentity(vertexGid[2]).id(),
                             vertexIndex);
    }

    Mpi_NewShadowSimplicesT() : simplexGid(empty)
    {
    }

    Mpi_NewShadowSimplicesT(Simplex *ref, int neiIndex, int myRank, unsigned int index) : simplexGid(empty)
    {
      set(ref, neiIndex, myRank, index);
    }

    bool isEmpty() const
    {
      return (simplexGid == empty);
    }

    void set(Simplex *ref, int neiIndex, int myRank, unsigned int index)
    {
      Simplex *nei = ref->getNeighbor(neiIndex);
      if (nei == NULL)
      {
        simplexGid = neighborIsVoid;
        __mpi_base_cell_index = index;
      }
      else
      {
        for (int i = 0; i < NVERT; i++) // NVERT==NNEI
        {
          vertexGid[i] = nei->getVertex(i)->getGlobalIdentity().get();
          if (nei->getNeighbor(i) == ref)
            vertexIndex = i;
        }
        Vertex *v = nei->getVertex(vertexIndex);
        generation = v->getGeneration().get();
        v->getCoords(vertexPos);
        simplexGid = nei->getGlobalIdentity(myRank).get();
        simplexGeneration = nei->getGeneration().get();
        __mpi_base_cell_index = index;
        vData = v->getData();
        sData = nei->getData();
      }
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<Coord>(+NDIM_W, OFFSETOF(MyType, vertexPos[0]));
      mpiDataType.push_back<GlobalIdentityValue>(+NVERT, OFFSETOF(MyType, vertexGid[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, simplexGid));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, simplexGeneration));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, generation));
      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, vertexIndex));
      mpiDataType.push_back<unsigned int>(1, OFFSETOF(MyType, __mpi_base_cell_index));

      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, sData));
      mpiDataType.push_back(VertexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, vData));

      mpiDataType.commit<MyType>();
      return mpiDataType;
    }
  };

  template <class T, class S>
  const typename Mpi_NewShadowSimplicesT<T, S>::GlobalIdentityValue Mpi_NewShadowSimplicesT<T, S>::empty =
      Mpi_NewShadowSimplicesT<T, S>::GlobalIdentity::empty.get();

  template <class T, class S>
  const typename Mpi_NewShadowSimplicesT<T, S>::GlobalIdentityValue Mpi_NewShadowSimplicesT<T, S>::neighborIsVoid =
      Mpi_NewShadowSimplicesT<T, S>::GlobalIdentity::max.get();

}

#include "../../internal/namespace.footer"
#endif
