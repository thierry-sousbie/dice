#ifndef __MPI_EXCHANGE_STRUCT_REPART_VERTICES__
#define __MPI_EXCHANGE_STRUCT_REPART_VERTICES__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  // unfortunately, Mpi_RefineQueryResultT has to be of type POD as it will be
  // used with MPI, so we cannot inherit and have to specialize the indexed
  // and the non indexed type for each and every MPI_ type ...
  template <class T, class V, bool WITH_INDEX = false>
  struct Mpi_RepartVerticesT
  {
    typedef Mpi_RepartVerticesT<T, V, WITH_INDEX> MyType;
    typedef V Vertex;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename T::GlobalIdentity::Value GlobalIdentityValue;
    typedef typename T::Coord Coord;
    typedef typename T::VertexData VertexData;
    static const int NDIM = T::NDIM;
    static const int NDIM_W = T::NDIM_W;
    static const GlobalIdentityValue empty;

    Coord coords[NDIM_W];
    GlobalIdentityValue generation;
    GlobalIdentityValue gid;
    VertexData vData;

    Mpi_RepartVerticesT()
    {
    }

    Mpi_RepartVerticesT(Vertex *v)
    {
      set(v);
    }

    bool isEmpty() const
    {
      return false;
    }

    void set(Vertex *v)
    {
      v->getCoords(coords);
      gid = v->getGlobalIdentity().get();
      generation = v->getGeneration().get();
      vData = v->getData();
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
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, generation));
      mpiDataType.push_back<GlobalIdentityValue>(1, OFFSETOF(MyType, gid));
      mpiDataType.push_back(VertexData::createDataElementsMpiStructType(), 1, OFFSETOF(MyType, vData));

      mpiDataType.commit<MyType>();
      return mpiDataType;
    }
  };

  template <class T, class V, bool WITH_INDEX>
  const typename Mpi_RepartVerticesT<T, V, WITH_INDEX>::GlobalIdentityValue Mpi_RepartVerticesT<T, V, WITH_INDEX>::empty =
      Mpi_RepartVerticesT<T, V, WITH_INDEX>::GlobalIdentity::empty.get();

} // namespace

#include "../../internal/namespace.footer"
#endif
