#ifndef __MPI_EXCHANGE_STRUCT_REFINE_SHARED__
#define __MPI_EXCHANGE_STRUCT_REFINE_SHARED__

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
  struct Mpi_RefineSharedT
  {
    typedef Mpi_RefineSharedT<T, WITH_INDEX> MyType;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef typename T::VertexData VertexData;
    typedef typename T::SimplexData SimplexData;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty; // = GlobalIdentity::empty.get();

    Coord coords[NDIM_W];                // coords of the new vertex
    GlobalIdentityValue segment[2];      // global id of the vertices of the split segment
    GlobalIdentityValue ref;             // global id of the reference simplex (the one that splits the segment)
    GlobalIdentityValue other;           // global ID of the partner simplex
    GlobalIdentityValue otherGeneration; // global ID of the partner simplex
    GlobalIdentityValue newVertexGeneration;
    GlobalIdentityValue newVertex; // global id of the newly created vertex
    SimplexData sData[2];
    VertexData vData;

    Mpi_RefineSharedT() : ref(empty)
    {
    }

    bool isEmpty() const
    {
      return (ref == empty);
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

      mpiDataType.push_back<Coord>(+NDIM_W, OFFSETOF(MyType, coords[0]));
      mpiDataType.push_back<GlobalIdentityValue>(2, OFFSETOF(MyType, segment[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, ref));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, other));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, otherGeneration));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, newVertexGeneration));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, newVertex));

      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 2, OFFSETOF(MyType, sData));
      mpiDataType.push_back(VertexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, vData));

      mpiDataType.commit<MyType>();

      return mpiDataType;
    }
  };

  template <class T>
  struct Mpi_RefineSharedT<T, true>
  {
    typedef Mpi_RefineSharedT<T, true> MyType;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef typename T::VertexData VertexData;
    typedef typename T::SimplexData SimplexData;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty; // = GlobalIdentity::empty.get();

    unsigned long __mpi_base_cell_index;
    Coord coords[NDIM_W];                // coords of the new vertex
    GlobalIdentity segment[2];           // global id of the vertices of the split segment
    GlobalIdentityValue ref;             // global id of the reference simplex (the one that split the segment)
    GlobalIdentityValue other;           // global ID of the partner simplex
    GlobalIdentityValue otherGeneration; // generation of the partner simplex
    GlobalIdentityValue newVertexGeneration;
    GlobalIdentityValue newVertex; // global id of the newly created vertex
    SimplexData sData[2];
    VertexData vData;

    Mpi_RefineSharedT() : ref(empty)
    {
    }

    bool isEmpty() const
    {
      return (ref == empty);
    }

    unsigned long getBaseCellIndex()
    {
      return __mpi_base_cell_index;
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<unsigned long>(1, OFFSETOF(MyType, __mpi_base_cell_index));
      mpiDataType.push_back<Coord>(+NDIM_W, OFFSETOF(MyType, coords[0]));
      mpiDataType.push_back<GlobalIdentityValue>(2, OFFSETOF(MyType, segment[0]));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, ref));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, other));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, otherGeneration));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, newVertexGeneration));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, newVertex));

      mpiDataType.push_back(SimplexData::createDataElementsMpiStructType(), 2, OFFSETOF(MyType, sData));
      mpiDataType.push_back(VertexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, vData));

      mpiDataType.commit<MyType>();

      return mpiDataType;
    }
  };

  template <class T, bool WITH_INDEX>
  const typename Mpi_RefineSharedT<T, WITH_INDEX>::GlobalIdentityValue Mpi_RefineSharedT<T, WITH_INDEX>::empty =
      Mpi_RefineSharedT<T, WITH_INDEX>::GlobalIdentity::empty.get();

  template <class T>
  const typename Mpi_RefineSharedT<T, true>::GlobalIdentityValue Mpi_RefineSharedT<T, true>::empty =
      Mpi_RefineSharedT<T, true>::GlobalIdentity::empty.get();

  template <class T, class S, bool WITH_INDEX = false>
  struct RefineSharedT : public Base<WITH_INDEX>
  {
    typedef Mpi_RefineSharedT<T, WITH_INDEX> MpiStruct;
    typedef S Simplex;
    typedef typename Simplex::Vertex Vertex;
    typedef typename Vertex::Data VertexData;

    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;

    Simplex *ref;
    Simplex *me;
    Simplex *other;
    Vertex *newVertex;
    Vertex *seg[2];
    int myRank;
    RefineSharedT() : ref(NULL)
    {
    }

    // DUMMY: this is not used for this type of structure
    void set(Simplex *curSimplex, const MpiStruct &str)
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("DUMMY function: RefineGhost::MpiStructs cannot be converted to RefineGhost directly.\n");
      glb::console->print<LOG_ERROR>("This function should never have been called.\n");
      exit(-1);
    }

    // ref_ is the simplex responsible for splitting the segment
    // me is a simplex that was split (i.e. the shrunk version of the former me)
    // other is the newly created partner of me (the other half of the former me)
    // newVertex is the new vertex introduced when splitting
    void set(Simplex *ref_, Simplex *me_, Simplex *other_, Vertex *newVertex_, int myRank_)
    {
      ref = ref_;
      me = me_;
      other = other_;
      newVertex = newVertex_;
      myRank = myRank_;

      // retrieve the 2 vertices that formed the segment which was split
      // -> they are the vertices that belong to 'me' or 'other' but not both
      for (int i = 0; i < Simplex::NVERT; i++)
      {
        if (other->getVertexIndex(me->getVertex(i)) < 0)
          seg[0] = me->getVertex(i);
        if (me->getVertexIndex(other->getVertex(i)) < 0)
          seg[1] = other->getVertex(i);
      }
    }

    void setEmpty()
    {
      ref = NULL;
    }

    bool isEmpty()
    {
      return (ref == NULL);
    }

    MpiStruct getMpiStruct()
    {
      MpiStruct result;

      if (ref != NULL)
      {
        Base<WITH_INDEX>::setMpiBaseCellIndex(result);
        result.ref = ref->getGlobalIdentity(myRank).get();
        result.other = other->getGlobalIdentity(myRank).get();
        result.otherGeneration = other->getGeneration().get();
        result.newVertex = newVertex->getGlobalIdentity().get();
        result.newVertexGeneration = newVertex->getGeneration().get();
        std::copy(newVertex->getCoordsPtr(),
                  newVertex->getCoordsPtr() + NDIM_W,
                  result.coords);
        result.segment[0] = seg[0]->getGlobalIdentity().get();
        result.segment[1] = seg[1]->getGlobalIdentity().get();

        result.sData[0] = me->getData();
        result.sData[1] = other->getData();

        result.vData = newVertex->getData(); //*static_cast<VertexData*>(newVertex);
      }

      return result;
    }
  };

} // namespace

#include "../../internal/namespace.footer"
#endif
