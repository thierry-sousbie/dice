#ifndef __MESH_HXX__
#define __MESH_HXX__

#include <set>

#include "../dice_globals.hxx"

#include "../mesh/localMesh.hxx"
#include "../mesh/meshParams.hxx"
#include "../mesh/mpiCellDataExchange.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_refineQueryResult.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_checkRefineResult.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_refineShared.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_repartSimplices.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_repartSimplicesWithCache.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_repartVertices.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_coarsenGhosts.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_coarsenShadows.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_simplicesGlobalIds.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_shadowSimplicesFromGhost.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_newShadowSimplices.hxx"

#include "../partition/partitioner.hxx"
#include "../partition/parmetisParams.hxx"

#include "../tools/types/cell.hxx"
#include "../tools/IO/myIO.hxx"
#include "../tools/IO/paramsParser.hxx"
#include "../tools/helpers/find_unordered_maps.hxx"
#include "../tools/sort/ompPSort.hxx"
#include "../tools/sort/peanoHilbert.hxx"

#include "../IO/ndNetUnstructuredMesh.hxx"

#include "../finiteElements/gaussQuadrature.hxx"

#include "internal/meshQuadratureFunctorAdapter.hxx"

/**
 * @file 
 * @brief  A class used to manage a MPI shared unstructured simplicial mesh 
 * embeded in a higher dimension space.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class MeshT
 * \brief A class used to manage a MPI shared unstructured simplicial mesh 
 * embeded in a higher dimension space. This class is responsible for managing 
 * global MPI operations, local operations being defined in the inherited class LocalMeshT
 * \tparam mesh traits. See MeshTraitsT.
 */

template <class T>
class MeshT : public LocalMeshT<T>
{
public:
  typedef MeshT<T> MyType;
  typedef LocalMeshT<T> LocalMesh;
  typedef typename LocalMesh::Tree Tree; 
  typedef typename Tree::AnyNodeBase TreeNode;

  static const int NDIM = LocalMesh::NDIM;
  static const int NDIM_W = LocalMesh::NDIM_W;  
  static const int BOUNDARY_TYPE = T::BOUNDARY_TYPE;
  static const int WORLD_BOUNDARY_TYPE = T::WORLD_BOUNDARY_TYPE; 
  static const int IS_PERIODIC = (BOUNDARY_TYPE==BoundaryType::PERIODIC);

  static std::string classHeader() {return "mesh";}
  static float classVersion() {return 0.11;}
  static float compatibleSinceClassVersion() {return 0.11;}

  typedef typename LocalMesh::Simplex Simplex;
  typedef typename LocalMesh::ShadowSimplex ShadowSimplex;
  typedef typename LocalMesh::GhostSimplex GhostSimplex;
  typedef typename LocalMesh::Facet Facet;
  typedef typename LocalMesh::Segment Segment;
  typedef typename LocalMesh::Vertex Vertex;
  typedef typename LocalMesh::ShadowVertex ShadowVertex;
  typedef typename LocalMesh::GhostVertex GhostVertex;

  typedef typename LocalMesh::SegmentHandle SegmentHandle; 
  typedef typename LocalMesh::FacetHandle FacetHandle;
 
  typedef typename T::Coord Coord;  
  typedef typename T::GlobalIndex    GlobalIndex;
  typedef typename T::LocalIndex     LocalIndex;
  typedef typename T::GlobalIdentity GlobalIdentity;
  
  typedef typename GlobalIdentity::Rank  GlobalIdentityRank;
  typedef typename GlobalIdentity::Id    GlobalIdentityId;
  typedef typename GlobalIdentity::Value GlobalIdentityValue;

  typedef typename LocalMesh::simplexPtr_iterator simplexPtr_iterator;
  typedef typename LocalMesh::ghostSimplexPtr_iterator ghostSimplexPtr_iterator;
  typedef typename LocalMesh::shadowSimplexPtr_iterator shadowSimplexPtr_iterator;
  typedef typename LocalMesh::vertexPtr_iterator  vertexPtr_iterator;
  typedef typename LocalMesh::ghostVertexPtr_iterator  ghostVertexPtr_iterator;
  typedef typename LocalMesh::shadowVertexPtr_iterator  shadowVertexPtr_iterator;

  typedef typename LocalMesh::const_simplexPtr_iterator const_simplexPtr_iterator;
  typedef typename LocalMesh::const_ghostSimplexPtr_iterator const_ghostSimplexPtr_iterator;
  typedef typename LocalMesh::const_shadowSimplexPtr_iterator const_shadowSimplexPtr_iterator;
  typedef typename LocalMesh::const_vertexPtr_iterator  const_vertexPtr_iterator;

  typedef typename LocalMesh::segment_circulator segment_circulator;
  typedef typename LocalMesh::const_segment_circulator const_segment_circulator;  

  typedef typename LocalMesh::GeometricProperties GeometricProperties;

  typedef typename LocalMesh::simplexPtr_LGS_iterator simplexPtr_LGS_iterator;
  typedef typename LocalMesh::vertexPtr_LGS_iterator vertexPtr_LGS_iterator;
  typedef typename LocalMesh::simplexPtr_LG_iterator simplexPtr_LG_iterator;
  typedef typename LocalMesh::vertexPtr_LG_iterator vertexPtr_LG_iterator;
  // typedef typename LocalMesh::const_allSimplexPtr_iterator const_allSimplexPtr_iterator;
  // typedef typename LocalMesh::const_allVertexPtr_iterator const_allVertexPtr_iterator;
  
  typedef typename Partitioner::Index PartitionerIndex;  

  // MPI structures
  // Used for coarsen/refine
  typedef mpiExchangeStruct::RefineQueryResultT<T,SegmentHandle,Simplex,true> 
  MpiExchg_RefineQueryResult;  
  typedef mpiExchangeStruct::RefineSharedT<T,Simplex,true> 
  MpiExchg_RefineShared;  
  typedef mpiExchangeStruct::Mpi_CoarsenGhostsT<T,Simplex> 
  MpiExchg_CoarsenGhosts; 
  typedef mpiExchangeStruct::Mpi_CoarsenShadowsT<T,Simplex> 
  MpiExchg_CoarsenShadows; 
  typedef mpiExchangeStruct::Mpi_NewShadowSimplicesT<T,Simplex> 
  MpiExchg_NewShadowSimplices;
  typedef mpiExchangeStruct::Mpi_CheckRefineResult 
  MpiExchg_CheckRefineResult;
  typedef mpiExchangeStruct::Mpi_SimplicesGlobalIdsT<T,Simplex> 
  MpiExchg_SimplicesGlobalIds;   

  // Used for repartitioning 
  typedef mpiExchangeStruct::Mpi_RepartSimplicesT<T,Simplex> 
  MpiExchg_RepartSimplices; 
  typedef mpiExchangeStruct::Mpi_RepartSimplicesWithCacheT<T,Simplex> 
  MpiExchg_RepartSimplicesWithCache;
  typedef mpiExchangeStruct::Mpi_RepartVerticesT<T,Vertex> 
  MpiExchg_RepartVertices;      
  typedef mpiExchangeStruct::Mpi_ShadowSimplicesFromGhostQueryT<T,Simplex> 
  MpiExchg_ShadowSimplicesFromGhostQuery;   
  typedef mpiExchangeStruct::Mpi_ShadowSimplicesFromGhostReplyT<T,Simplex> 
  MpiExchg_ShadowSimplicesFromGhostReply;   

  typedef CellT<> GridCell;
  typedef SimplicialGridT< NDIM, NDIM_W,GridCell,BOUNDARY_TYPE==BoundaryType::PERIODIC> 
  SimplicialGrid;
  //typedef ImplicitTesselationT<NDIM,NDIM_W,GridCell> ImplicitTesselation;  
  typedef MeshParamsT<NDIM,NDIM_W,Coord> Params; 
  
  typedef typename LocalMesh::CheckRefineReturnType CheckRefineReturnType;

  MeshT(MpiCommunication *com=glb::mpiComWorld)
    :LocalMesh(com),
     mpiCom(com),
     ghostExchange(com),
     shadowExchange(com) 
  {
    LocalMesh::construct();
    construct();
  }
     

private:

  void construct()
  {
    // Allocate mpi types used for data transfer
    mpiType_refineQueryResult=
      MpiExchg_RefineQueryResult::MpiStruct::createMpiStructType();    
    mpiType_refineShared=
      MpiExchg_RefineShared::MpiStruct::createMpiStructType();
    mpiType_checkRefineResult=
      MpiExchg_CheckRefineResult::createMpiStructType();
    mpiType_coarsenGhosts=
      MpiExchg_CoarsenGhosts::createMpiStructType();
    mpiType_coarsenShadows=
      MpiExchg_CoarsenShadows::createMpiStructType();
    mpiType_newShadowSimplices=
      MpiExchg_NewShadowSimplices::createMpiStructType();
    mpiType_simplicesGlobalIds=
      MpiExchg_SimplicesGlobalIds::createMpiStructType();

    // used for repartitionning
    mpiType_repartSimplices=
      MpiExchg_RepartSimplices::createMpiStructType();
    mpiType_repartSimplicesWithCache=
      MpiExchg_RepartSimplicesWithCache::createMpiStructType();
    mpiType_repartVertices=
      MpiExchg_RepartVertices::createMpiStructType();
    mpiType_shadowSimplicesFromGhostQuery=
      MpiExchg_ShadowSimplicesFromGhostQuery::createMpiStructType();
    mpiType_shadowSimplicesFromGhostReply=
      MpiExchg_ShadowSimplicesFromGhostReply::createMpiStructType();
   
    //std::copy(&params->x0[0],&params->x0[NDIM_W],std::back_inserter(x0));
    //std::copy(&params->delta[0],&params->delta[NDIM_W],std::back_inserter(delta));
    #pragma omp critical
    {   
      mpiTagsStart_repart = mpiCom->reserveTags(10);
      mpiTagsStart_refine = mpiCom->reserveTags(5);
      mpiTagsStart_coarsen = mpiCom->reserveTags(10);
      mpiTagsStart_general = mpiCom->reserveTags(5);
    }
  }
  
public:

  ~MeshT()
  {  
    //MPI_Type_free(&mpiType_refineQueryResult);
    //MPI_Type_free(&mpiType_checkRefineResult);
    //MPI_Type_free(&mpiType_coarsenShadows);
    //MPI_Type_free(&mpiType_simplicesGlobalIds);
    
    //delete params;
    mpiCom->barrier();
  }
  
  /** \brief Defragment the mesh memory pools and restrict memory usage to the strict 
   *  minimum. 
   *
   * Calling this function may be usefull after a large number of elements has been 
   * refined and then coarsened. Note that defragmenting may take time as the mesh must 
   * be temporarily written to disk. 
   */
  void defrag()
  {
    LocalMesh::defrag();
  }

  /** \brief serialize the mesh to a file. 
   *  \param writer a pointer to the writer
   *  \tparam W a writer class such as myIO::BinaryWriterT
   */
  template <class W>
  void write(W *writer)
  {    
    writer->writeHeader(classHeader(),classVersion());
    params.write(writer);    
    
    LocalMesh::write(writer);
    
    writer->write(globalNLocalCells,NDIM+1);
    writer->write(globalNCellsCum,NDIM+1);
    writer->write(&loadImbalanceFactor);    

    ghostExchange.serialize(writer);
    shadowExchange.serialize(writer);    
  }

  /** \brief unserialize the mesh from a file written using write().
   *  \param writer a pointer to the writer
   *  \tparam W a writer class such as myIO::BinaryWriterT
   */
  template <class R>
  void read(Params &meshParams, R *reader)
  {    
    //if (reader == NULL) return build(meshParams);
    
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);

    params = meshParams;
    params.read(reader);

    if (glb::console->willPrint<LOG_INFO>())
      {
	glb::console->print<LOG_STD>("Unserializing the mesh:\n");
	glb::console->indent();
      }
    else
      glb::console->printFlush<LOG_STD>("Unserializing the mesh ... ");

    typename LocalMesh::UVPUpdater vpu;
    typename LocalMesh::UGVPUpdater gvpu;
    typename LocalMesh::USVPUpdater svpu;

    typename LocalMesh::USPUpdater spu;
    typename LocalMesh::UGSPUpdater gspu;
    typename LocalMesh::USSPUpdater sspu;   

    LocalMesh::build(params,reader,vpu,gvpu,svpu,spu,gspu,sspu);
 
    reader->read(globalNLocalCells,NDIM+1);
    reader->read(globalNCellsCum,NDIM+1);
    reader->read(&loadImbalanceFactor);
  
    ghostExchange.unSerialize(reader, *spu, *gspu);
    shadowExchange.unSerialize(reader, *spu, *sspu);

    if (glb::console->willPrint<LOG_INFO>())
      {
	glb::console->unIndent();
	glb::console->print<LOG_STD>("Done.\n");	
      }
    else
    glb::console->print<LOG_STD>("done.\n");

    if (glb::debug) 
      {
	checkConsistencyAndReport<LOG_ERROR>("restart");
	ghostExchange.template print<LOG_STD_ALL>("ghosts");
	shadowExchange.template print<LOG_STD_ALL>("shadow");	
      }

    updateCellsCount();
    glb::console->printToBuffer<LOG_INFO>("The tesselation has %ld vertices and %ld simplices.\n",getGlobalNCells(0),getGlobalNCells(NDIM));
    glb::console->flushBuffer<LOG_INFO>();

    
    //mpiCom->barrier();
    //exit(0);
  }

  /** \brief Builds the mesh from an implicit tesselation.
   *  \param implicitTesselation the desrired tesselation
   *  \param meshParams the mesh parameters
   *  \tparam IST A class describing an implicit tesselation conforming to the interface
   *  defined by ImplicitTesselationT
   */
  template <class IST>
  void build(const IST *implicitTesselation, const Params &meshParams)
  { 
    params = meshParams;
    LocalMesh::build(meshParams);

    glb::console->printNewLine<LOG_STD>
      ("Creating initial mesh from implicit tesselation:\n");
    glb::console->indent();
    /*
    SimplicialGrid sg(params.x0,params.delta,
		      params.resolution,
		      params.initTesselationType);
    */
    std::vector<unsigned long> nCells=implicitTesselation->getNCells();
    /*
    glb::console->printToBuffer<LOG_INFO>
      ("Initial grid resolution is [%ld",(long)params.resolution[0]);
    for (int i=1;i<NDIM;i++) 
      glb::console->printToBuffer<LOG_INFO>(",%ld",(long)params.resolution[i]);
    glb::console->printToBuffer<LOG_INFO>("].\n");
    glb::console->flushBuffer<LOG_INFO>();
    */

    glb::console->printToBuffer<LOG_INFO>("BBox dimensions: x0=[%ld",(long)params.x0[0]);
    for (int i=1;i<NDIM;i++) 
      glb::console->printToBuffer<LOG_INFO>(",%ld",(long)params.x0[i]);
    glb::console->printToBuffer<LOG_INFO>("], delta=[%ld",(long)params.delta[0]);
    for (int i=1;i<NDIM;i++) 
      glb::console->printToBuffer<LOG_INFO>(",%ld",(long)params.delta[i]);
    glb::console->printToBuffer<LOG_INFO>("].\n");
    glb::console->flushBuffer<LOG_INFO>();

    glb::console->printToBuffer<LOG_INFO>("Tesselation has %ld vertices and [%ld",
					  nCells[0],nCells[1]);
    for (int i=2;i<=NDIM;i++) 
      glb::console->printToBuffer<LOG_INFO>(",%ld",nCells[i]);
    glb::console->printToBuffer<LOG_INFO>("] k-cells.\n");
    glb::console->flushBuffer<LOG_INFO>();

    std::vector< PartitionerIndex > initialPartition;
    if (mpiCom->size()>1)
      {
	glb::console->printFlush<LOG_STD>("Building %d partitions ... ",mpiCom->size());
	long ncut=Partitioner::buildFromGrid(implicitTesselation,
					     initialPartition,params.initPartitionType,
					     params.initPartitionTolerance);
	if (ncut>=0) glb::console->print<LOG_STD>("done. (%ld cuts)\n",ncut);
	else glb::console->print<LOG_STD>("done.\n");
      }
    
    //glb::memoryInspector->report<LOG_INFO>();    

    glb::console->print<LOG_STD>("Initializing local %ldD mesh (%ldD space):\n",
				 NDIM,NDIM_W);
    glb::console->indent();

    initFromImplicitTesselation(implicitTesselation,initialPartition);
    initialPartition.clear();    

    glb::console->unIndent();   
    glb::console->print<LOG_STD>("DONE.\n");
    
    glb::memoryInspector->reportAllProcesses<LOG_INFO>();

    glb::console->unIndent();
    glb::console->print<LOG_STD>("Creation DONE.\n");
   
    // char fname[80];
    // sprintf(fname,"test_%4.4d.NDnet",mpiCom->rank());
    // dumpToNDnetwork(fname);       
  }
  
  /** \brief Save the mesh to a NDnet file. 
   *  \param fname the name of the file
   *  \param options options given as ::IO::NDNET_WriterOptions
   *  \param completeFName set to true to add a number corresponding to the rank of the 
   *  MPI process at the end of fname, so that each individual pieces have a different 
   *  filename.
   */
  void dumpToNDnetwork(const std::string &fname, 
		       int options=IO::NDNET_Default,
		       bool completeFName=false)
  {
    typedef IO::NDnetUnstructuredMeshWriterT<MyType> NdNetWriter;
    std::string name(fname);
    
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d",fname.c_str(),mpiCom->rank());
	name = std::string(tmp);
      }  
    /*
    else
      {
	char tmp[255];
	sprintf(tmp,"%s",fname.c_str());
	name = std::string(tmp);
      }
    */

    // long options=IO::NDNET_Default;
    // if (withNeighbors) options |= IO::NDNET_WithNeighbors;
    // if (dumpAll) options |= IO::NDNET_WithShadows | IO::NDNET_WithGhosts;

    NdNetWriter ndNetWriter(this,options,name.c_str());
    ndNetWriter.write();
    mpiCom->barrier();
  }

  /** \brief Save the mesh to a NDnet file after filtering it.
   * \param fname the name of the file
   * \param filter1 a functor taking a \a CT1* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter2 a functor taking a \a CT2* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter3 a functor taking a \a CT3* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter4 a functor taking a \a CT4* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param options options given as ::IO::NDNET_WriterOptions
   * \param completeFName set to true to add a number corresponding to the rank of the 
   * MPI process at the end of fname, so that each individual pieces have a different 
   * filename.
   * \tparam F1 The type of the 1st filter functor
   * \tparam F2 The type of the 2nd filter functor
   * \tparam F3 The type of the 3rd filter functor
   * \tparam F4 The type of the 4th filter functor
   * \tparam CT1 The type of cell filtered by the 1st filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT2 The type of cell filtered by the 2nd filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT3 The type of cell filtered by the 3rd filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT4 The type of cell filtered by the 4th filter (Simplex, FacetHandle or SegmentHandle)
   */
  template <class F1, class F2, class F3, class F4,
	    class CT1=typename F1::CellType,
	    class CT2=typename F2::CellType,
	    class CT3=typename F3::CellType,
	    class CT4=typename F4::CellType>
  void dumpToFilteredNDNetwork(const std::string &fname,
			       const F1 &filter1,
			       const F2 &filter2,
			       const F3 &filter3,
			       const F4 &filter4,
			       int options=IO::NDNET_Default,
			       bool completeFName=false)
  {
    typedef IO::NDnetUnstructuredMeshWriterT<MyType> NdNetWriter;
    std::string name(fname);
    
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d",fname.c_str(),mpiCom->rank());
	name = std::string(tmp);
      }  
    
    NdNetWriter ndNetWriter(this,options,name.c_str());
    ndNetWriter.template filter<F1,CT1>(filter1);
    ndNetWriter.template filter<F2,CT2>(filter2);
    ndNetWriter.template filter<F3,CT3>(filter3);
    ndNetWriter.template filter<F4,CT4>(filter4);
    ndNetWriter.write();
  }

  /** \brief Save the mesh to a NDnet file after filtering it.
   * \param fname the name of the file
   * \param filter1 a functor taking a \a CT1* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter2 a functor taking a \a CT2* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter3 a functor taking a \a CT3* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.   
   * \param options options given as ::IO::NDNET_WriterOptions
   * \param completeFName set to true to add a number corresponding to the rank of the 
   * MPI process at the end of fname, so that each individual pieces have a different 
   * filename.
   * \tparam F1 The type of the 1st filter functor
   * \tparam F2 The type of the 2nd filter functor
   * \tparam F3 The type of the 3rd filter functor  
   * \tparam CT1 The type of cell filtered by the 1st filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT2 The type of cell filtered by the 2nd filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT3 The type of cell filtered by the 3rd filter (Simplex, FacetHandle or SegmentHandle)  
   */
  template <class F1, class F2, class F3,
	    class CT1=typename F1::CellType,
	    class CT2=typename F2::CellType,
	    class CT3=typename F3::CellType>
  void dumpToFilteredNDNetwork(const std::string &fname,
			       const F1 &filter1,
			       const F2 &filter2,
			       const F3 &filter3,
			       int options=IO::NDNET_Default,
			       bool completeFName=false)
  {
    typedef IO::NDnetUnstructuredMeshWriterT<MyType> NdNetWriter;
    std::string name(fname);
    
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d",fname.c_str(),mpiCom->rank());
	name = std::string(tmp);
      }  
    
    NdNetWriter ndNetWriter(this,options,name.c_str());
    ndNetWriter.template filter<F1,CT1>(filter1);
    ndNetWriter.template filter<F2,CT2>(filter2);
    ndNetWriter.template filter<F3,CT3>(filter3);
    
    ndNetWriter.write();
  }

  /** \brief Save the mesh to a NDnet file after filtering it.
   * \param fname the name of the file
   * \param filter1 a functor taking a \a CT1* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.
   * \param filter2 a functor taking a \a CT2* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.   
   * \param options options given as ::IO::NDNET_WriterOptions
   * \param completeFName set to true to add a number corresponding to the rank of the 
   * MPI process at the end of fname, so that each individual pieces have a different 
   * filename.
   * \tparam F1 The type of the 1st filter functor
   * \tparam F2 The type of the 2nd filter functor   
   * \tparam CT1 The type of cell filtered by the 1st filter (Simplex, FacetHandle or SegmentHandle)
   * \tparam CT2 The type of cell filtered by the 2nd filter (Simplex, FacetHandle or SegmentHandle)   
   */
  template <class F1, class F2,
	    class CT1=typename F1::CellType,
	    class CT2=typename F2::CellType>
  void dumpToFilteredNDNetwork(const std::string &fname,
			       const F1 &filter1,
			       const F2 &filter2,
			       int options=IO::NDNET_Default,
			       bool completeFName=false)
  {
    typedef IO::NDnetUnstructuredMeshWriterT<MyType> NdNetWriter;
    std::string name(fname);
    
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d",fname.c_str(),mpiCom->rank());
	name = std::string(tmp);
      }  
    
    NdNetWriter ndNetWriter(this,options,name.c_str());
    ndNetWriter.template filter<F1,CT1>(filter1);
    ndNetWriter.template filter<F2,CT2>(filter2);    
    ndNetWriter.write();
  }
  
  /** \brief Save the mesh to a NDnet file after filtering it.
   * \param fname the name of the file
   * \param filter1 a functor taking a \a CT1* object \a cell and returning NULL or \a cell 
   * to exclude / include it in the output file.   
   * \param options options given as ::IO::NDNET_WriterOptions
   * \param completeFName set to true to add a number corresponding to the rank of the 
   * MPI process at the end of fname, so that each individual pieces have a different 
   * filename.
   * \tparam F1 The type of the 1st filter functor 
   * \tparam CT1 The type of cell filtered by the 1st filter (Simplex, FacetHandle or SegmentHandle)     
   */
  template <class F,class CT=typename F::CellType>
  void dumpToFilteredNDNetwork(const std::string &fname,
			       const F &filter,
			       int options=IO::NDNET_Default,
			       bool completeFName=false)
  {
    typedef IO::NDnetUnstructuredMeshWriterT<MyType> NdNetWriter;
    std::string name(fname);
    
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d",fname.c_str(),mpiCom->rank());
	name = std::string(tmp);
      }  
    
    NdNetWriter ndNetWriter(this,options,name.c_str());
    ndNetWriter.template filter<F,CT>(filter);
    ndNetWriter.write();
  }

  // std::vector<unsigned long> getGlobalNCells()
  // {
  //   return globalNCellsCum[i].back()
  //   return globalNCells;
  // }

  /** \brief Get cells count in the global mesh
   *  \param type the type of cell
   *  \return the number of cells of type \a type in the global mesh
   */
  unsigned long getGlobalNCells(int type)
  {
    return (globalNCellsCum[type].size()==0)?0:globalNCellsCum[type].back();
    //return globalNCells[i];
  }

  /** \brief Synchronize the value stored in the cache variable of ghost simplices to 
   *  the value stored in the cache variable of their local simplex.   
   */
  template <typename CacheT = uint64_t>
  void synchronizeGhostSimplicesCache()
  {        
    std::vector< std::vector<CacheT> > send(ghostExchange.sendRank.size());
    std::vector< std::vector<CacheT> > receive;

#pragma omp parallel for
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];
	send[i].resize(ghostExchange.send[src].size());
	for (unsigned long j=0;j<send[i].size();j++)
	  memcpy(&send[i][j],&(ghostExchange.send[src][j]->cache.ui64),sizeof(CacheT));
	//send[i][j] = ghostExchange.send[src][j]->cache.ui64;
      }
    
    ghostExchange.exchange(send,receive,0);
    
#pragma omp parallel for
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int src = ghostExchange.receiveRank[i];
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    GhostSimplex *g=ghostExchange.receive[src][j];
	    memcpy(&(g->cache.ui64),&receive[i][j],sizeof(CacheT));
	    //g->cache.ui64 = receive[i][j];
	  }
      }			   
  }
  
  /** \brief Synchronize the value stored in the cache variable of shadow simplices to 
   *  the value stored in the cache variable of their local simplex.   
   */
  template <typename CacheT = uint64_t>
  void synchronizeShadowSimplicesCache()
  {
    std::vector< std::vector<CacheT> > send(shadowExchange.sendRank.size());
    std::vector< std::vector<CacheT> > receive;

#pragma omp parallel for
    for (unsigned long i=0;i<shadowExchange.sendRank.size();i++)
      {
	int src = shadowExchange.sendRank[i];
	send.resize(shadowExchange.send[src].size());
	for (unsigned long j=0;j<send[i].size();j++)
	  memcpy(&send[i][j],&(shadowExchange.send[src][j]->cache.ui64),sizeof(CacheT));
	//send[i][j] = shadowExchange.send[src][j]->cache.ui64;
      }
    
    shadowExchange.exchange(send,receive,0);

#pragma omp parallel for
    for (unsigned long i=0;i<shadowExchange.receiveRank.size();i++)
      {
	int src = shadowExchange.receiveRank[i];
	for (unsigned long j=0;j<shadowExchange.receive[src].size();j++)
	  {
	    ShadowSimplex *g=shadowExchange.receive[src][j];
	    memcpy(&(g->cache.ui64),&receive[i][j],sizeof(CacheT));
	    //g->cache.ui64 = receive[i][j];
	  }
      }			   
  }
 
  /** \brief Sort the simplices and vertices of the mesh locally (i.e. independantly on 
   *  each MPI process) according to an order defined by the \a sf and \a vf functors.
   *  \param sf A functor taking a Simplex* as argument and returning an object of a type 
   *  comparable with operator<(..) to the one returned by \a vf.
   *  \param vf A functor taking a Vertex* as argument and returning an object of a type 
   *  comparable with operator<(..) to the one returned by \a sf.
   *  \param nThreads The number of openMP threads to use for parallel sort
   */
  template <class SimplexFunctor, class VertexFunctor>
  void sort(const SimplexFunctor &sf, 
	    const VertexFunctor &vf,
	    int nThreads=glb::num_omp_threads)
  {    
    typename LocalMesh::SVPUpdater vpu;
    typename LocalMesh::SGVPUpdater gvpu;
    typename LocalMesh::SSVPUpdater svpu;

    typename LocalMesh::SSPUpdater spu;
    typename LocalMesh::SGSPUpdater gspu;
    typename LocalMesh::SSSPUpdater sspu; 

    LocalMesh::sort(sf,vf,spu,gspu,sspu,vpu,gvpu,svpu,nThreads);
    
    ghostExchange.updatePointers(spu,gspu);
    shadowExchange.updatePointers(spu,sspu); 
    // Simplices pointers have changed, so incidence became invalid !
    incidentSimplices.needFullUpdate=true;
  }

  /*
  template <class ET>
  GlobalIndex getGlobalIndex(ET *element)
  {
    return globalIdentityToGlobalIndex(element->getGlobalIdentity(),element->getType());
  }
  */

  /** \brief Repartition the mesh in order to improve load balance if needed.
   *
   *  The imbalance is computed as  \f$ i=max_k(weight[k]*nSimplices[k]) / 
   *  avg_k(weight[k]*nSimplices[k]) \f$ where \a k is the rank of each process. 
   *  Repartitionning only happens of i>params.repartThreshold (see MeshParamsT).
   *  \param weight A relative weight given to the local partition before
   *  computing the load balance. If \a weight is 0, all partitions have the same weight.
   *  \param force if true, forces repartitionning to happen
   */
  // FIXME : post an Irecv before Isend and use waitall ...
  // FIXME : it would be nice NOT to reallocate a new ghostSimplex pool ...
  // FIXME : can we replace mpi_all2all by something more local ??
  // FIXME : add a function to clean the Queue in memoryPool ?
  // Weight is the weight of this process, if weight<=0, then the weight
  // is given by the number of local cells
  bool repart(double weight=0, bool force=false, int nThreads=glb::num_omp_threads)
  {    
    typedef typename my_dense_set<Vertex*>::type VertexDenseSet;
    typedef typename VertexDenseSet::iterator VertexDenseSet_it; 

    typedef typename my_dense_set<Simplex*>::type SimplexDenseSet;
    typedef typename SimplexDenseSet::iterator SimplexDenseSet_it; 

    typedef typename my_dense_set<GlobalIdentityValue>::type GidDenseSet;
    typedef typename GidDenseSet::iterator GidDenseSet_it; 

    typedef typename my_dense_hash<GlobalIdentityValue,Vertex*>::type VertexDenseHash;
    typedef typename VertexDenseHash::iterator VertexDenseHash_it; 

    typedef typename my_dense_hash<GlobalIdentityValue,Simplex*>::type SimplexDenseHash;
    typedef typename SimplexDenseHash::iterator SimplexDenseHash_it;         

    if (mpiCom->size()==1) return false;

    double imbalance;
    double weightPerCell;

    if (weight<=0) 
      {
	weightPerCell=1.0;
	imbalance = getLoadImbalanceFactor();
      }
    else
      {
	weightPerCell = weight / LocalMesh::getNCells()[NDIM];

	auto minMaxSum = mpiCom->minMaxSum(weight);
	double max=minMaxSum.first.second;
	double avg=minMaxSum.second / mpiCom->size();
	imbalance = max/avg;
      }

    if ((imbalance < params.repartThreshold)&&(!force)) 
      {
	glb::console->printFlush<LOG_INFO>("Repartitioning skipped (f=%.3f <= %.2f).\n",
					   imbalance,
					   params.repartThreshold);
	return false;  
      }

    // WE NEED TO REPARTITION !
  
    // -> incidence vectors are now invalid !
    clearIncidentSimplices();
    
    if (glb::debug>1) 
      dumpToNDnetwork("before_repart",
		      IO::NDNET_WithShadows | 
		      IO::NDNET_WithGhosts,
		      true); 
      //dumpToNDnetwork("before_repart",true,true,false); 
   
    const int myRank = mpiCom->rank();
    const int nParts = mpiCom->size();
    std::vector< std::vector<MPI_Request> > requests(6);    
    
    if (glb::console->willPrint<LOG_INFO>())
      glb::console->print<LOG_INFO>("Repartitioning (f=%.3f):\n",imbalance);
    else
      glb::console->printFlush<LOG_STD>("Repartitioning (f=%.3f) ... ",imbalance);
    glb::console->indent();

    std::vector<PartitionerIndex> partition;
    MpiCellDataExchangeT<Simplex,TreeNode> leavesExchange(mpiCom);
    Tree::repart(partition,
		 params.refinePartitionType,
		 params.repartTolerance,
		 weightPerCell,
		 leavesExchange,
		 nThreads);

    // FIXME : show this ?
    //leavesExchange.template print<LOG_DEBUG>("leaves");   
    leavesExchange.template print<LOG_PEDANTIC>("leaves");

    glb::console->printFlush<LOG_INFO>("Rebuilding the mesh ... ");

    const unsigned long nVertices=LocalMesh::vertexPool.getUsedCount();
    //const unsigned long nSimplices=LocalMesh::simplexPool.getUsedCount();

    // this is used to allocate vectors to store simplices
    // that contain a given vertex
    typedef MemoryPoolT< std::vector<Simplex*> , true> SimplexVectorPool;
    SimplexVectorPool simplexVectorPool("vector<Simplex*>",2.0);
    simplexVectorPool.reserve(100000);

    // vertexBallArr stores simplices containing a given vertex.
    // Only vertices on interface will be used so this is OK memory wise
    // FIXME: use fixed size elements pool instead of vector ?
    std::vector<std::vector<Simplex*> *> vertexBallArr(nVertices,NULL);
  
    // Set the simplices cache to the rank of the process 
    // they will belong to after repartitionning
    LocalMesh::setSimplicesCacheL(myRank); // default -> remain local
    
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<leavesExchange.sendRank.size();++i)
      {
	const long index = leavesExchange.sendRank[i];
	const std::vector<Simplex *> &curSend = leavesExchange.send[index];	
	for (unsigned long j=0;j<curSend.size();++j)
	  curSend[j]->cache.l=index;
      }

    // receiveRankIndex is used to retrieve the index in leaveExchange corresponding
    // to a given process rank
    std::vector<int> receiveRankIndex(nParts,-1);
    for (unsigned long i=0;i<leavesExchange.receiveRank.size();++i)
      receiveRankIndex[leavesExchange.receiveRank[i]]=i;

    // We also need to retrieve the destination of each ghost simplices
    // and store it in the cache
    std::vector< std::vector<int> > ghostSimplexDestSnd(ghostExchange.sendRank.size());
    //std::vector<MPI_Request> ghostRequests(ghostExchange.sendRank.size());
    requests[0].resize(ghostExchange.sendRank.size());
    // send the data ...   
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
      {
	const std::vector<Simplex *> &curSend = 
	  ghostExchange.send[ghostExchange.sendRank[i]];

	ghostSimplexDestSnd[i].resize(curSend.size());
	for (unsigned long j=0;j<curSend.size();++j)
	  ghostSimplexDestSnd[i][j] = curSend[j]->cache.l;
	//#pragma omp critical
	mpiCom->Isend(&ghostSimplexDestSnd[i][0],ghostSimplexDestSnd[i].size(),
		      ghostExchange.sendRank[i],&requests[0][i],mpiTagsStart_repart+0);
      }
   
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
      {	
	int source=mpiCom->Probe(mpiTagsStart_repart+0); 
	std::vector<int> received(ghostExchange.receive[source].size());
	//received.resize(ghostExchange.receive[source].size());
	mpiCom->Recv(&received[0],received.size(),source,mpiTagsStart_repart+0);
	for (unsigned long j=0;j<received.size();++j)
	  ghostExchange.receive[source][j]->cache.l = received[j];
      }
   
    // for (int i=0;i<ghostExchange.sendRank.size();++i)
    //   mpiCom->wait(&ghostRequests[i]);
    // ghostSimplexDestSnd.clear();
        
    // ISend and IReceive the simplices that need to be exchanged in the background,
    // so that we can cover the latency with some computations
    // we also build sendVertexSet, to store which vertices will need to be sent
    
    // first the Ireceive part    
    std::vector<MPI_Request> receivedSimplicesRequests(leavesExchange.receiveRank.size());
    std::vector< std::vector<MpiExchg_RepartSimplices> > 
      receivedSimplices(leavesExchange.receiveRank.size());
    for (unsigned long i=0;i<leavesExchange.receiveRank.size();++i)
      {	
	int source=leavesExchange.receiveRank[i];
	receivedSimplices[i].resize(leavesExchange.receive[source].size());	
	mpiCom->IrecvMpiType(&receivedSimplices[i][0],
			     receivedSimplices[i].size(),
			     mpiType_repartSimplices.getType(),source,
			     &receivedSimplicesRequests[i],mpiTagsStart_repart+1);
      }   
  
    // And then the Isend part 
    std::vector< std::vector<MpiExchg_RepartSimplices> > 
      sendSimplices(leavesExchange.sendRank.size()); 
    // +1 => last vector stores the vertices kept locally !
    std::vector< VertexDenseSet > sendVertexSet(leavesExchange.sendRank.size()+1);
    set_set_empty_key(sendVertexSet.back(),NULL);
    // local hash
    const long sendVertexSetBackLoadGuess = 0.2*LocalMesh::getNVertices();
    sendVertexSet.back().rehash(sendVertexSetBackLoadGuess/
				(sendVertexSet.back().max_load_factor()-0.01)); 
    long sendVertexSetLoadGuess[leavesExchange.sendRank.size()];
    requests[1].resize(leavesExchange.sendRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<leavesExchange.sendRank.size();++i)
      {
	const std::vector<Simplex *> &curSend = 
	  leavesExchange.send[leavesExchange.sendRank[i]];
	
	sendSimplices[i].resize(curSend.size());
	set_set_empty_key(sendVertexSet[i],NULL);
	//set_set_deleted_key(sendVertexSet[i],GlobalIdentity::max.get());

	sendVertexSetLoadGuess[i]=1.5*curSend.size();
	sendVertexSet[i].rehash(sendVertexSetLoadGuess[i]/ 
				(sendVertexSet[i].max_load_factor()-0.01));
	for (unsigned long j=0;j<curSend.size();++j)
	  {
	    sendSimplices[i][j].set(curSend[j],myRank);
	    for (int k=0;k<Simplex::NVERT;k++)	      
	      sendVertexSet[i].insert(curSend[j]->getVertex(k));
	  }
	//#pragma omp critical
	mpiCom->IsendMpiType(&sendSimplices[i][0],sendSimplices[i].size(),
			     mpiType_repartSimplices.getType(),
			     leavesExchange.sendRank[i],&requests[1][i],
			     mpiTagsStart_repart+1);
      }
    
    // Identify the newly shared vertices, that will be at the interface 
    // of different processes after repartionning.
    // To do that, we identify neighbor simplices with different destinations.      
    //FIXME : this should be openmp friendly with a critical section on the pop() ...
    //#pragma omp parallel for
    //for (unsigned long i=0;i<glb::num_omp_threads;i++)
    //{
    const simplexPtr_iterator its_end=LocalMesh::simplexEnd();
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)
      {

	for (int k=0;k<Simplex::NNEI;++k)
	  {
	    Simplex *nei=it->getNeighbor(k);
	    if (nei==NULL) continue;
	    
	    // if (it->cache.l != nei->cache.l), simplices go to different places
	    // => the common face will be an interface between different processes
	    if (it->cache.l != nei->cache.l)
	      {	
		for (int j=0;j<Simplex::NVERT;++j)
		  {
		    if (j==k) continue; // this one is NOT on the interface
		    
		    long lid=it->getVertex(j)->getLocalIndex();
		    if (vertexBallArr[lid]==NULL)
		      {
			// pop a vector to store the simplices surrounding this vertex
#pragma omp critical
			simplexVectorPool.pop(&vertexBallArr[lid]);
			//FIXME: estimate the expected number of simplices ?
			// vertexBallArr[lid]->reserve((1<<NDIM)*2); 
		      }		  
		  }
	      }	   
	  }
	
      }
    
    // Don't forget to scan the ghosts as we may miss some of the local shared 
    // vertices if we don't
    const ghostSimplexPtr_iterator itgs_end=LocalMesh::ghostSimplexEnd();
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)    
      {

	for (int k=0;k<Simplex::NNEI;++k)
	  {
	    Simplex *nei=it->getNeighbor(k);
	    if (nei==NULL) continue;
	    
	    // if (it->cache.l != nei->cache.l), simplices go to different places
	    // => the common face will be an interface between different processes
	    if (it->cache.l != nei->cache.l)
	      {	
		for (int j=0;j<Simplex::NVERT;++j)
		  {
		    if (it->getVertex(j)->isShadowOrGhost()) continue;
		    if (j==k) continue; // this one is NOT on the interface
		    
		    long lid=it->getVertex(j)->getLocalIndex();
		    if (vertexBallArr[lid]==NULL)
		      {
			// pop a vector to store the simplices surrounding this vertex
#pragma omp critical
			simplexVectorPool.pop(&vertexBallArr[lid]);
			//FIXME: estimate the expected number of simplices ?
			// vertexBallArr[lid]->reserve((1<<NDIM)*2); 
		      }		  
		  }
	      }	   
	  }
	
      }

    // We now want to retrieve the simplices surrounding interface vertices   
    // first allocate some space to store them    
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)	 
      {
	for (int j=0;j<Simplex::NVERT;++j)
	  {	  
	    long lid=it->getVertex(j)->getLocalIndex();
	    if (vertexBallArr[lid]!=NULL)
	      vertexBallArr[lid]->push_back(*it);
	  }	
      }

    // don't forget the ghosts ...
    //ghostSimplexPtr_iterator itgs_end=LocalMesh::ghostSimplexEnd();
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)	 
      {
	for (int j=0;j<Simplex::NVERT;++j)
	  {
	    Vertex *v=it->getVertex(j);
	    if (v->isShadowOrGhost()) continue;
	    long lid=v->getLocalIndex();
	    
	    if (vertexBallArr[lid]!=NULL)
	      vertexBallArr[lid]->push_back(*it);
	  }
      }

    // We can now identify the new ghosts simplices for the local process
    // as well as for the processes to which we will send simplices
    std::vector<int> sendRankIndex(nParts,-1);    
    std::vector<SimplexDenseSet> newGhostSimplices(leavesExchange.sendRank.size()+1);
    long newGhostSimplicesLoadGuess[leavesExchange.sendRank.size()];
    for (unsigned long i=0;i<leavesExchange.sendRank.size();++i)
      {
	set_set_empty_key(newGhostSimplices[i],NULL);
	sendRankIndex[leavesExchange.sendRank[i]]=i;
	newGhostSimplicesLoadGuess[i]=3*sendSimplices[i].size();
	newGhostSimplices[i].rehash(newGhostSimplicesLoadGuess[i]/
				    (newGhostSimplices[i].max_load_factor()-0.01));
      }
    set_set_empty_key(newGhostSimplices.back(),NULL);  
    const long newGhostSimplicesBackLoadGuess = 0.1*LocalMesh::getNSimplices();
    newGhostSimplices.back().rehash(newGhostSimplicesBackLoadGuess/
				    (newGhostSimplices.back().max_load_factor()-0.01)); 
  

    // Each interface vertex has at least 2 simplices with different
    // destinations. Let P be the rank of the destination process of one of the
    // simplices surrounding an interface vertex V, then any simplex containing V
    // with a destination rank Q != P belongs the the ghost layer of P.
    // simplex->cache.l contains the destination rank (new local rank after repart)
    vertexPtr_iterator itv_end=LocalMesh::vertexEnd();
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
      {       
	long lid=it->getLocalIndex();
	// we will need clean tags soon, so do that now ...
	it->cleanTags(); 
	if (vertexBallArr[lid] != NULL)
	  {
	    // first retrieve the set of processes ranks the simplices
	    // surrounding this vertex will go to
	    std::set<long> dest;
	    for (unsigned long i=0;i<vertexBallArr[lid]->size();++i)
	      {
		Simplex* cur=vertexBallArr[lid]->at(i);
		// we don't want ghosts as they will be taken care of on the 
		// process they belong ...
		if (!cur->isShadowOrGhost()) dest.insert(cur->cache.l);	
	      }

	    // For each rank in the set, all the simplices with different new
	    // rank will be ghost
	    for (typename std::set<long>::iterator dst_it=dest.begin();
		 dst_it!=dest.end();++dst_it)
	      {	
		// N.B.: new ghosts of the local process are stored in the
		// back of newGhostsSimplices
		int idx=(*dst_it == myRank)?
		  (newGhostSimplices.size()-1):
		  sendRankIndex[*dst_it];		

		for (unsigned long i=0;i<vertexBallArr[lid]->size();++i)
		  {
		    Simplex *s=vertexBallArr[lid]->at(i);
		    if (s->cache.l == *dst_it) continue; 
		    
		    newGhostSimplices[idx].insert(s);
		    for (int j=0;j<Simplex::NVERT;j++)
		      {		
			sendVertexSet[idx].insert(s->getVertex(j));
			//Vertex *v=s->getVertex(j);
			//sendVertexHash[idx].insert(std::make_pair(v->getGlobalIdentity().get(),v));
		      }
		  }
	      }	      
	  }
      }

     // send the vertices that need to be sent
    std::vector< std::vector<MpiExchg_RepartVertices> > 
      sendVertices(leavesExchange.sendRank.size());
    requests[3].resize(leavesExchange.sendRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<leavesExchange.sendRank.size();++i)
      {	
	sendVertices[i].resize(sendVertexSet[i].size());
	unsigned long j=0;
	const VertexDenseSet_it it_end=sendVertexSet[i].end();
	for (VertexDenseSet_it it=sendVertexSet[i].begin();it!=it_end;++it)
	  sendVertices[i][j++].set(*it);

	// glb::console->print<LOG_DEBUG>("sendVertexSet[%ld]: size = %ld  => %.2g \n",i,
	// 			       sendVertexSet[i].size(),
	// 			       double(sendVertexSet[i].size())/
	// 			       sendVertexSetLoadGuess[i]);
	sendVertexSet[i].clear();

	mpiCom->IsendMpiType(&sendVertices[i][0],sendVertices[i].size(),
			     mpiType_repartVertices.getType(),
			     leavesExchange.sendRank[i],&requests[3][i],
			     mpiTagsStart_repart+3);
      }    
    
    // and send the new remote ghosts. Meanwhile, we can do some work. 
    requests[2].resize(leavesExchange.sendRank.size());
    std::vector< std::vector<MpiExchg_RepartSimplicesWithCache> > 
      sendGhostSimplices(leavesExchange.sendRank.size());    

#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<leavesExchange.sendRank.size();++i)
      {
	//mpiCom->wait(&requests[1][i]); // make sure we can clear sendSimplices
	sendGhostSimplices[i].resize(newGhostSimplices[i].size());
	unsigned long j=0;
	const SimplexDenseSet_it it_end=newGhostSimplices[i].end();
	for (SimplexDenseSet_it it=newGhostSimplices[i].begin();it!=it_end;++it)
	  sendGhostSimplices[i][j++].set(*it,myRank);
	
	// glb::console->print<LOG_DEBUG>("newGhostSimplices[%ld]: size = %ld  => %.2g \n",i,
	// 			       newGhostSimplices[i].size(),
	// 			       double(newGhostSimplices[i].size())/
	// 			       newGhostSimplicesLoadGuess[i]);
	newGhostSimplices[i].clear(); //don't need that anymore !

	mpiCom->IsendMpiType(&sendGhostSimplices[i][0],sendGhostSimplices[i].size(),
			     mpiType_repartSimplicesWithCache.getType(),
			     leavesExchange.sendRank[i],&requests[2][i],
			     mpiTagsStart_repart+2);
      }
         
    // the new ghosts that are already local need to be serialized as 
    // if they had been sent via MPI ...
    long counter=0;
    std::vector<MpiExchg_RepartSimplicesWithCache> 
      newLocalGhostSimplices(newGhostSimplices.back().size());

    for (SimplexDenseSet_it it=newGhostSimplices.back().begin();
	 it != newGhostSimplices.back().end();++it)
      newLocalGhostSimplices[counter++].set(*it,myRank);
    

    glb::console->print<LOG_DEBUG>("newGhostSimplicesBack: size = %ld  => %.2g \n",
				   newGhostSimplices.back().size(),
				   double(newGhostSimplices.back().size())/
				   newGhostSimplicesBackLoadGuess);
    newGhostSimplices.clear(); //don't need that anymore !
    
    // start building a hash table with simplices and vertices that will 
    // remain but whose neighborhood/vertices will require updating
    // Here we also tag the vertices on the interface of the simplices
    // that remain local and add them to the hash table
    const long simplicesHashLoadGuess = 0.25*LocalMesh::getNSimplices();
    const long verticesHashLoadGuess = 0.25*LocalMesh::getNVertices();
    const long localSimplicesReserveGuess = 0.03*LocalMesh::getNSimplices();

    SimplexDenseHash simplicesHash;
    set_hash_empty_key(simplicesHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(simplicesHash,GlobalIdentity::max.get());
    simplicesHash.rehash(simplicesHashLoadGuess/
			 (simplicesHash.max_load_factor()-0.01));
    VertexDenseHash verticesHash;
    set_hash_empty_key(verticesHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(verticesHash,GlobalIdentity::max.get());
    verticesHash.rehash(verticesHashLoadGuess/
			(verticesHash.max_load_factor()-0.01));
    std::vector<MpiExchg_RepartSimplices> localSimplices;
    localSimplices.reserve(localSimplicesReserveGuess);
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)	 
      {
	if (it->cache.l != myRank) //this simplex won't stay local
	  continue; // note that we can't recycle it yet, as we need the neighbors !
	/*
	// we will store the boundary layer in the hash later, using the tagged2
	// vertices, so no need to add it here ...
	// simplicesHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	bool boundaryLayer=false;
	bool pushed=false;
	for (int k=0;k<Simplex::NVERT;++k)
	  {
	    Vertex *v=it->getVertex(k);
	    long lid=v->getLocalIndex();
	    // Select only those vertices that will become shared but local
	    if (vertexBallArr[lid]!=NULL)
	      {
		boundaryLayer=true;
		if (!v->isTagged2())
		  {
		    v->setTagged2F();
		    verticesHash.insert(std::make_pair(v->getGlobalIdentity().get(),v));
		  }
	      }
	  }
	if (boundaryLayer)
	  {
	    MpiExchg_RepartSimplices mpiSimplex(*it,myRank);
	    localSimplices.push_back(mpiSimplex);
	    pushed=true;
	  }
	for (int j=0;j<Simplex::NNEI;++j)
	  {
	    Simplex *nei=it->getNeighbor(j);
	    if (nei==NULL) continue;
	    
	    // check if the neighbor is/will be on a different process
	    if (nei->isShadowOrGhost())
	      {	     
		if (!pushed)
		  {
		    MpiExchg_RepartSimplices mpiSimplex(*it,myRank);
		    localSimplices.push_back(mpiSimplex);
		    pushed=true;
		  }
		for (int k=0;k<Simplex::NVERT;++k)
		  {
		    if (j==k) continue; // this one is NOT on the interface !
		    Vertex *v=it->getVertex(k);
		    // tag2 the interface vertices to allow recovering the local 
		    // simplices boundary layer 
		    if (!v->isTagged2())
		      {
			v->setTagged2F();
			verticesHash.insert(std::make_pair(v->getGlobalIdentity().get(),v));
		      }
		  }
	      }		     
	  }
	*/
	  
	for (int j=0;j<Simplex::NNEI;++j)
	  {
	    Simplex *nei=it->getNeighbor(j);
	    if (nei==NULL) continue;
	    
	    // check if the neighbor is/will be on a different process
	    if ((it->cache.l != nei->cache.l)||(nei->isShadowOrGhost()))
	      {	     
		MpiExchg_RepartSimplices mpiSimplex(*it,myRank);
		localSimplices.push_back(mpiSimplex);
		// we will store the boundary layer in the hash later, using the tagged2
		// vertices, so no need to add it here ...
		// simplicesHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
		for (int k=0;k<Simplex::NVERT;++k)
		  {
		    if (j==k) continue; // this one is NOT on the interface !
		    Vertex *v=it->getVertex(k);
		    // tag2 the interface vertices to allow recovering the local 
		    // simplices boundary layer 
		    if (!v->isTagged2())
		      {
			v->setTagged2F();
			verticesHash.insert(std::make_pair(v->getGlobalIdentity().get(),v));
		      }
		  }
	      }		     
	  }
	  	  
      }

    // We must not forget to recycle allocated vertexBalls before destroying the pool
    for (unsigned long i=0;i<vertexBallArr.size();++i)
      if (vertexBallArr[i]!=NULL)
	simplexVectorPool.recycle(&vertexBallArr[i]);
    vertexBallArr.clear();
    
    // Recycle simplices that are NOT local anymore so that we can use
    // them to store the new local simplices, and tag the vertices of 
    // simplices that remain local.
    // Interface vertices were tagged2 so we can also add to the simplexHash
    // all the simplices on the inner boundary layer of the set of local
    // simplices (needed to reconnect the transfered pieces)    
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)	 
      {
	if (it->cache.l != myRank)
	  LocalMesh::simplexPool.recycle(*it);
	else
	  {
	    for (int j=0;j<Simplex::NVERT;++j)
	      {
		Vertex *v=it->getVertex(j);
		v->setTaggedF(); 
		if (v->isTagged2()) // interface vertex
		  simplicesHash.insert(std::make_pair(it->getGlobalIdentity(myRank).get(),
						      *it));
	      }
	  }	  
      }

    // time to receive the new simplices that were sent earlier
    // we do it with Isend/Ireceive earlier
    mpiCom->Waitall(receivedSimplicesRequests);
    /*
    // This is the version where we don't use Irecv
    std::vector<int> receiveRankIndex(nParts,-1);
    for (long i=0;i<leavesExchange.receiveRank.size();++i)
    receiveRankIndex[leavesExchange.receiveRank[i]]=i;

    std::vector< std::vector<MpiExchg_RepartSimplices> > 
       receivedSimplices(leavesExchange.receiveRank.size());
    for (int i=0;i<leavesExchange.receiveRank.size();++i)
    {	
    int source=mpiCom->Probe(1); 
    int srcIndex=receiveRankIndex[source];
    receivedSimplices[srcIndex].resize(leavesExchange.receive[source].size());	
    mpiCom->RecvMpiType(&receivedSimplices[srcIndex][0],
    receivedSimplices[srcIndex].size(),
    mpiType_repartSimplices,source,1);
    }    
    */
    
    // allocate the newly received simplices and add them to the hash table
    for (unsigned long i=0;i<leavesExchange.receiveRank.size();++i)
      {
	long index = leavesExchange.receiveRank[i];
	for (unsigned long j=0;j<leavesExchange.receive[index].size();++j)
	  { 
	    Simplex *simplex;
	    LocalMesh::simplexPool.pop(&simplex);	
	    simplex->setData(receivedSimplices[i][j].sData);
	    simplex->setGeneration(receivedSimplices[i][j].generation);
	    leavesExchange.receive[index][j]->addChild(simplex);	     
	    simplicesHash.insert(std::make_pair(receivedSimplices[i][j].gid,simplex));
	  }
      }  
 
    // we will allocate the new ghosts in a temporary pool
    // that we can swap with the current one.
    // As a bonus, the ghost pool will be defragmented ...
    typename LocalMesh::GhostSimplexPool newGhostSimplexPool("new GhostSimplex");
    newGhostSimplexPool.reserve(LocalMesh::ghostSimplexPool.getUsedCount());
  
    // Allocate and add the local ghosts to the hash
    for (unsigned long i=0;i<newLocalGhostSimplices.size();++i)
      {
	GhostSimplex *gSimplex;
	newGhostSimplexPool.pop(&gSimplex);
	gSimplex->setData(newLocalGhostSimplices[i].sData);
	gSimplex->setGeneration(newLocalGhostSimplices[i].generation);
	simplicesHash.insert(std::make_pair(newLocalGhostSimplices[i].gid,
					    gSimplex));
      }

    // Now we can receive the new non-local ghosts.
    // We then allocate and add them to the hash only if they are not there already.      
    std::vector< std::vector<MpiExchg_RepartSimplicesWithCache> > 
      receivedGhostSimplices(leavesExchange.receiveRank.size());

    for (unsigned long i=0;i<leavesExchange.receiveRank.size();++i)
      {	
	int count=0;
	int source=mpiCom->ProbeCount(count,mpiType_repartSimplicesWithCache.getType(),
				      mpiTagsStart_repart+2);
	int srcIndex=receiveRankIndex[source];

	if (count>0) receivedGhostSimplices[srcIndex].resize(count);
	else receivedGhostSimplices[srcIndex].clear();
	  
	mpiCom->RecvMpiType(&receivedGhostSimplices[srcIndex][0],count,
			    mpiType_repartSimplicesWithCache.getType(),source,
			    mpiTagsStart_repart+2);
	
	for (int j=0;j<count;++j)
	  {	    
	    GlobalIdentityValue gid=receivedGhostSimplices[srcIndex][j].gid;
	    const SimplexDenseHash_it it=simplicesHash.find(gid);
	    if (it==simplicesHash.end())
	      {
		GhostSimplex *gSimplex;
		newGhostSimplexPool.pop(&gSimplex);
		gSimplex->setData(receivedGhostSimplices[srcIndex][j].sData);
		gSimplex->setGeneration(receivedGhostSimplices[srcIndex][j].generation);
		simplicesHash.insert(std::make_pair(gid,gSimplex));
	      }
	  }	
      }  
   
    // the already local vertices will be either rebuilt if they are ghosts/shadow
    // or conserved if not
    std::vector<MpiExchg_RepartVertices> newLocalVertices;
    newLocalVertices.reserve(sendVertexSet.back().size());
    for (VertexDenseSet_it it=sendVertexSet.back().begin();
	 it!=sendVertexSet.back().end();++it)
      {
	Vertex *v=*it;
	if (v->isShadowOrGhost()) 
	  {
	    // if the vertex is shadow/ghost, it needs to be poped again in 
	    // a different pool, so we store its data before recycling it
	    //MpiExchg_RepartVertices tmp(v);
	    newLocalVertices.push_back(MpiExchg_RepartVertices(v));	    	 
	  }
	else
	  {
	    // the vertex is already in the vertex pool
	    v->setTaggedF();
	    verticesHash.insert(std::make_pair(v->getGlobalIdentity().get(),v));
	  }
      }	

    glb::console->print<LOG_DEBUG>("sendVertexSetBack: size = %ld  => %.2g \n",
				   sendVertexSet.back().size(),
				   double(sendVertexSet.back().size())/
				   sendVertexSetBackLoadGuess);
    sendVertexSet.clear();

    // Now we can recycle the non tagged vertices (they are not local anymore) 
    // and remove any previous Shared flag
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
      {
	if (!it->isTagged()) {LocalMesh::vertexPool.recycle(*it);}
	else it->setSharedF(false);
      }
    
    //glb::console->print<LOG_PEDANTIC_ALL>("after : %ld vertices.\n",LocalMesh::getNVertices());
   
    // Here only the local vertices that remained local are in the pool, and their shared
    // flags have been wiped -> receive the new vertices !
    for (unsigned long i=0;i<leavesExchange.receiveRank.size();++i)
      {	
	int count=0;
	int source=mpiCom->ProbeCount(count,mpiType_repartVertices.getType(),
				      mpiTagsStart_repart+3);	
	//int srcIndex=receiveRankIndex[source];
	std::vector<MpiExchg_RepartVertices> receivedVertices(count);		
	mpiCom->RecvMpiType(&receivedVertices[0],
			    receivedVertices.size(),
			    mpiType_repartVertices.getType(),source,
			    mpiTagsStart_repart+3);
	for (unsigned long j=0;j<receivedVertices.size();++j)
	  {   
	    MpiExchg_RepartVertices &mpiVertex=receivedVertices[j];
	    if (verticesHash.find(mpiVertex.gid)==verticesHash.end())
	      {
		Vertex *vertex;
		LocalMesh::vertexPool.pop(&vertex);
		verticesHash.insert(std::make_pair(mpiVertex.gid,vertex));
		vertex->setCoords(mpiVertex.coords);	
		vertex->setData(mpiVertex.vData);	
		vertex->setGlobalIdentity(mpiVertex.gid);
		vertex->setGeneration(mpiVertex.generation);
		vertex->setSetF();	
	      }
	  }
      } 

    // and create the vertices we are missing !
     for (unsigned long i=0;i<newLocalVertices.size();++i)
      {	
	if (verticesHash.find(newLocalVertices[i].gid)==verticesHash.end())
	  {
	    Vertex *vertex;
	    LocalMesh::vertexPool.pop(&vertex);
	    verticesHash.insert(std::make_pair(newLocalVertices[i].gid,vertex));
	    
	    vertex->setCoords(newLocalVertices[i].coords);
	    vertex->setGlobalIdentity(newLocalVertices[i].gid);
	    vertex->setGeneration(newLocalVertices[i].generation);
	    vertex->setData(newLocalVertices[i].vData);
	    vertex->setSetF();
	  }
      }  

    //glb::console->print<LOG_PEDANTIC_ALL>("after2 : %ld vertices.\n",LocalMesh::getNVertices());

    // Now we can rebuild everything, going through the simplices in :
    //    - newLocalGhostSimplices : the already local new ghosts
    //    - receivedGhostSimplices : the non-local new ghosts
    //    - receivedSimplices : the simplices imported from another process
    //    - localSimplices: the local simplices whose neighbors changed  
    // If a neighbor global id is not found from the hash, this means it is
    // a shadow simplex
    std::vector< std::vector<MpiExchg_RepartSimplices> * > allNewSimplices;
    for (unsigned long i=0;i<receivedSimplices.size();++i) 
      allNewSimplices.push_back(&receivedSimplices[i]); 
    allNewSimplices.push_back(&localSimplices); // This must be the last one !
    
    // first set the new local simplices
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<allNewSimplices.size();++i)
      {
	const std::vector<MpiExchg_RepartSimplices> &cur = *allNewSimplices[i];
	for (unsigned long j=0;j<cur.size();++j)
	  {	
	    Simplex *s=simplicesHash.find(cur[j].gid)->second;
	  
	    // set the vertices if that's not done already 
	    if (!s->isSet())
	      {
		Vertex *v[Simplex::NVERT];
		for (int k=0;k<Simplex::NVERT;++k)
		  v[k]=verticesHash.find(cur[j].vertices[k])->second;
		s->setVertices(v);
		s->setSafeF();
	      }

	    for (int k=0;k<Simplex::NNEI;++k)
	      {
		if (cur[j].neighbors[k] == MpiExchg_RepartSimplices::empty)
		  {
		    s->neighbors[k]=NULL;
		    s->setBoundaryF();
		    for (int l=0;l<Simplex::NVERT;++l)
		      if (l!=k) 
			{
			  #pragma omp critical
			  s->getVertex(l)->setBoundaryF();
			}
		  }
		else
		  {
		    const SimplexDenseHash_it it=simplicesHash.find(cur[j].neighbors[k]);
		    if (it!=simplicesHash.end()) 
		      s->neighbors[k]=it->second;
		    //else if (i!=allNewSimplices.size()-1) // the neighbor is a shadow
		  }
	      }	    
	  }
      }

    // and then the ghost simplices
    std::vector< std::vector<MpiExchg_RepartSimplicesWithCache> * > allNewGhostSimplices;
    allNewGhostSimplices.push_back(&newLocalGhostSimplices);
    for (unsigned long i=0;i<receivedGhostSimplices.size();++i) 
      allNewGhostSimplices.push_back(&receivedGhostSimplices[i]);
    
    // we will store the shadows in a hash so that we can retrieve them
    SimplexDenseHash shadowSimplicesHash;
    set_hash_empty_key(shadowSimplicesHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(shadowSimplicesHash,GlobalIdentity::max.get());
    const long shadowSimplicesHashLoadGuess = LocalMesh::getNShadowSimplices();
    shadowSimplicesHash.rehash(shadowSimplicesHashLoadGuess/
			       (shadowSimplicesHash.max_load_factor()-0.01));

    // because all the shadow simplices will be reallocated, we free the previous ones
    const shadowSimplexPtr_iterator itss_end=LocalMesh::shadowSimplexEnd();
    for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();it!=itss_end;++it)
      LocalMesh::shadowSimplexPool.recycle(*it);        
    
    //FIXME: not so easy to make openMP parallel, as the same simplex may be represented 
    // more than once ...
    //FIXME #pragma omp parallel for
    for (unsigned long i=0;i<allNewGhostSimplices.size();++i)
      {
	const std::vector<MpiExchg_RepartSimplicesWithCache> &cur=*allNewGhostSimplices[i];
	for (unsigned long j=0;j<cur.size();++j)
	  {	
	    Simplex *s=simplicesHash.find(cur[j].gid)->second;
	    if (s->isGhost()) // should be always true 
	      {
		static_cast<GhostSimplex *>(s)->setGlobalIdentity(cur[j].gid);
		static_cast<GhostSimplex *>(s)->setGeneration(cur[j].generation);
		s->cache.l = cur[j].cache; // that's the original copy's new rank
	      }

	    // set the vertices if that's not done already 
	    if (!s->isSet())
	      {
		Vertex *v[Simplex::NVERT];
		for (int k=0;k<Simplex::NVERT;++k)
		  {
		    //VertexDenseHash_it it=verticesHash.find(cur[j].vertices[k]);
		    v[k]=verticesHash.find(cur[j].vertices[k])->second;
		  }
		s->setVertices(v);
		s->setSafeF();
	      }

	    for (int k=0;k<Simplex::NNEI;++k)
	      {
		if (cur[j].neighbors[k] == MpiExchg_RepartSimplices::empty)
		  {
		    s->neighbors[k]=NULL;
		    s->setBoundaryF();
		    if (!s->isShadow())
		      {
			for (int l=0;l<Simplex::NVERT;++l)
			  if (l!=k) s->getVertex(l)->setBoundaryF();
		      }
		  }
		else
		  {
		    const SimplexDenseHash_it it=simplicesHash.find(cur[j].neighbors[k]);
		    if (it!=simplicesHash.end()) 
		      s->neighbors[k]=it->second;
		    else // this is a shadow
		      {	
#pragma omp critical
			{ // START OF CRITICAL 			  		  
			  SimplexDenseHash_it it=
			    shadowSimplicesHash.find(cur[j].neighbors[k]);
			  ShadowSimplex *shadow;
			  if (it==shadowSimplicesHash.end())
			    {			     
			      // shadow neighbor does not exist, create it.
			      LocalMesh::shadowSimplexPool.pop(&shadow);
			      shadow->setLocalIndex(LocalMesh::getNShadowSimplices()-1);
			      //shadow->setGlobalIdentity(cur[j].neighbors[k]);	
			      shadowSimplicesHash.insert
				(std::make_pair(cur[j].neighbors[k],
						static_cast<Simplex*>(shadow)));
			    }
			  else shadow=static_cast<ShadowSimplex*>(it->second);
			  s->neighbors[k]=shadow;
			} // END OF CRITICAL 
					
		      }	//else (it!=simplicesHash.end()) 	
		  } //else (cur[j].neighbors[k] == MpiExchg_RepartSimplices::empty)
	      }	//for k
	  }//for j
      }//for i

    glb::console->print<LOG_DEBUG>("shadowSimplicesHash: size = %ld  => %.2g \n",
				   shadowSimplicesHash.size(),
				   double(shadowSimplicesHash.size())/
				   shadowSimplicesHashLoadGuess);
    shadowSimplicesHash.clear();

    // we now have to swap the ghostSimplexPool with the new one so that we can 
    // iterate over the new ghost simplices. Note that we need to recycle all the 
    // old ghosts as their destructor won't be called when the pool get destroyed.    
    //glb::console->print<LOG_STD_ALL>("before :%ld\n",LocalMesh::getNGhostSimplices());
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
      LocalMesh::ghostSimplexPool.recycle(*it);
    LocalMesh::ghostSimplexPool.swap(newGhostSimplexPool);
    
    // clean all the tags from local vertices 
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)	 
      it->cleanTags(); 

    // Now some ghosts vertices are actually allocated from simplexPool rather than
    // ghostSimplexPool, so we need to correct that and set shared vertices flags
    // so we tag any vertex belonging to a ghost simplex
    //const ghostSimplexPtr_iterator itgs_end=LocalMesh::ghostSimplexEnd();    
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();
	 it!=LocalMesh::ghostSimplexEnd();++it)
      {
	for (int i=0;i<Simplex::NVERT;++i)
	  it->getVertex(i)->setTagged2F();
      }
    
    // a vertex is shared if it belongs to a ghost and a local simplex at least
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)	 
      {
	for (int i=0;i<Simplex::NVERT;++i)
	  if (it->getVertex(i)->isTagged2())
	    {	      
	      it->getVertex(i)->cleanTags(); // so that only the ghost vertices are tagged2
	      it->getVertex(i)->setSharedF();
	    }
      }   
   
    // We can now identify ghost vertices and reallocate them
    VertexDenseHash ghostVertexHash;
    set_hash_empty_key(ghostVertexHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(ghostVertexHash,GlobalIdentity::max.get());
    const long ghostVertexHashLoadGuess = LocalMesh::getNGhostVertices();
    ghostVertexHash.rehash(ghostVertexHashLoadGuess/
			   (ghostVertexHash.max_load_factor()-0.01));

    // recycle all previous ghost vertices
    const ghostVertexPtr_iterator itgv_end=LocalMesh::ghostVertexEnd();
    for (ghostVertexPtr_iterator it=LocalMesh::ghostVertexBegin();it!=itgv_end;++it)
      {
	LocalMesh::ghostVertexPool.recycle(*it);
      }

    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
      {
	int gvCount=0;
	for (int j=0;j<Simplex::NVERT;++j)
	  {
	    Vertex *v=it->getVertex(j);
	    
	    // shared vertices have been untagged !
	    if (v->isTagged2())
	      {
		VertexDenseHash_it ghost_it=
		  ghostVertexHash.find(v->getGlobalIdentity().get());
		GhostVertex *gv;
		if (ghost_it==ghostVertexHash.end())
		  {		    
		    LocalMesh::ghostVertexPool.pop(&gv);
		    gv->copy(v); // this will preserve the ghost status (flags) of gv
		    ghostVertexHash.insert(std::make_pair(gv->getGlobalIdentity().get(),
							  gv));
		  }
		else gv = static_cast<GhostVertex*>(ghost_it->second);
		it->setVertex(j,gv);
		gv->setGlobalIdentity(GlobalIdentity::max);
		gvCount++;
		// ghost vertex new global identity will be retrieved later !
		// as the lowest ranked global ID of its copies, so we ensure
		// it will have a higher rank than any copy.
	      }
	  }
	
	// Check that ghost simplices having more than a single shared vertex do not have
	// shadow simplex neighbors as this is not possible !
	if (gvCount<NDIM)
	  {
	    bool failed=false;
	    for (int i=0;i<Simplex::NNEI;++i) 
	      if(it->getNeighbor(i)->isShadow()) failed=true;
	    if (failed)
	      {
		glb::console->print<LOG_STD_ALL>
		  ("WARNING: This ghost simplex may not have a shadow neighbor !\n");
		static_cast<Simplex*>(*it)->
		  template print<LOG_STD_ALL>("Ghost:");
		for (int n=0;n<Simplex::NVERT;++n)
		  (*it)->getVertex(n)->template print<LOG_STD_ALL>();
	      }
	  }
	
      }
  
    // remove the new ghost vertices from the vertexPool and clean-up the tags
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
      {
	if (it->isTagged2())
	  LocalMesh::vertexPool.recycle(*it);
	it->cleanTags();
      }

    // Now the local IDs should change, as they need to be strictly growing from 0
    // so we update them locally ... 
    long lid=0;
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
      it->setLocalIndex(lid++);

    lid=0;
    for (ghostVertexPtr_iterator it=LocalMesh::ghostVertexBegin();it!=itgv_end;++it)
      {
	it->setLocalIndex(lid++);
	it->cleanTags();
      }

    // Note that this will also update the global IDs of the simplices !
    lid=0;
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)
      it->setLocalIndex(lid++);
    
    // Now we need to retrieve the orginals of ghost simplices
    std::vector<int> ghostOrigin(nParts,0);        

    lid=0;
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
      {
	it->setLocalIndex(lid++);
	ghostOrigin[it->cache.l]++; // count how many ghosts come from which rank
      }

    // we also need to rebuild the ghostExchange and shadowExchange structures
    // but for now we only know where the GhostSimplices originals are stored
    ghostExchange.receiveRank.clear();    
    for (unsigned long i=0;i<ghostOrigin.size();i++)
      {
	ghostExchange.receive[i].clear();
	if (ghostOrigin[i]>0)
	  {
	    ghostExchange.receiveRank.push_back(i);
	    ghostExchange.receive[i].reserve(ghostOrigin[i]);	   
	  }
      }
  
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
      ghostExchange.receive[it->cache.l].push_back(*it);

    // Now we want to send the ghosts global IDs so that we can retrieve 
    // their new one and identify them to their "original" (i.e. the one
    // version that is local on its process)
    // First we need to know how many requests each process will receive    
    std::fill(ghostOrigin.begin(),ghostOrigin.end(),0);
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
      {
	int src=ghostExchange.receiveRank[i];
	// ghostOrigin will contain how many simplices we request from each rank
	ghostOrigin[src]=ghostExchange.receive[src].size(); 
      }
    
    int requestsCount = 0;
    // FIXME : use allreduce or alltoall ????
    // mpiCom->Allreduce_inplace(ghostOrigin,MPI_SUM);
    // requestsCount = ghostOrigin[myRank];    
    std::vector<int> ghostSendStats(nParts);
    // after Alltoall we know how many simplices each rank wants from us
    mpiCom->Alltoall(ghostOrigin,ghostSendStats);
    for (unsigned long i=0;i<ghostSendStats.size();++i)
      if (ghostSendStats[i]>0) requestsCount++;
    
    //glb::console->print<LOG_STD_ALL>("Requests count : %d\n",requestsCount);

    // Send the global Ids of the requested ghosts
    std::vector< std::vector<GlobalIdentityValue> > 
      ghostGid(ghostExchange.receiveRank.size());

    requests[4].resize(ghostExchange.receiveRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
      {
	int src=ghostExchange.receiveRank[i];
	const std::vector<GhostSimplex *> &curSend = ghostExchange.receive[src];
	ghostGid[i].resize(curSend.size());
	for (unsigned long j=0;j<curSend.size();j++)
	  ghostGid[i][j]=curSend[j]->getGlobalIdentity().get();
	mpiCom->Isend(&ghostGid[i][0],ghostGid[i].size(),
		      ghostExchange.receiveRank[i],&requests[4][i],mpiTagsStart_repart+4);
      }
  
    // now we can allocate the ghostExchange.send
    ghostExchange.sendRank.clear();    
    for (unsigned long i=0;i<ghostSendStats.size();i++)
      {
	ghostExchange.send[i].clear();
	if (ghostSendStats[i]>0)
	  {
	    ghostExchange.sendRank.push_back(i);
	    ghostExchange.send[i].reserve(ghostSendStats[i]);	   
	  }
      }

    // and receive the global ids of the ghosts simplices we have to send
    for (int i=0;i<requestsCount;i++)
      {
	int source=mpiCom->Probe(mpiTagsStart_repart+4); 
	std::vector<GlobalIdentityValue> received(ghostSendStats[source]);
	mpiCom->Recv(&received[0],received.size(),source,mpiTagsStart_repart+4);
	for (unsigned long j=0;j<received.size();++j)
	  {
	    SimplexDenseHash_it it= simplicesHash.find(received[j]);
	    if (it == simplicesHash.end())
	      {
		PRINT_SRC_INFO(LOG_ERROR);
		glb::console->print<LOG_ERROR>
		  ("Could not retrieve ghost simplex origin.\n");
		exit(-1);
	      }
	    ghostExchange.send[source].push_back(it->second);
	  }	
      }
 
    ghostExchange.updateNCum();   
    // We can now update the local vertices global IDs (this is already done for simplices)
    for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
      it->setGlobalIdentity(myRank,it->getLocalIndex());      
   
    std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > 
      sendGlobalIds(ghostExchange.sendRank.size());
    std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > receivedGlobalIds;
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];
	sendGlobalIds[i].resize(ghostExchange.send[src].size());
	for (unsigned long j=0;j<sendGlobalIds[i].size();j++)
	  sendGlobalIds[i][j].set(ghostExchange.send[src][j],myRank);
      }

    ghostExchange.exchangeStruct(sendGlobalIds,receivedGlobalIds,
				 mpiType_simplicesGlobalIds.getType(),true,
				 mpiTagsStart_repart+6);
    
    // and set the new global Ids ...    
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int src = ghostExchange.receiveRank[i];
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    GhostSimplex *g=ghostExchange.receive[src][j];
	    MpiExchg_SimplicesGlobalIds &mpiData=receivedGlobalIds[i][j];
	    for (int k=0;k<Simplex::NVERT;k++)
	      {
		Vertex *v=g->getVertex(k);
		GlobalIdentity globalIdentity(mpiData.vertices[k]);
		// local shared vertices (i.e. those on the interface) 
		// have the global ID of the copy with lowest rank
		if (v->getGlobalIdentity().rank()>globalIdentity.rank())
		  v->setGlobalIdentity(globalIdentity);
	      }
	    g->setGlobalIdentity(mpiData.gid);
	  }
      }
    
    // Now we need to transmit again the shared vertices, as when they also are 
    // ghost vertices of a remote node, the remote node does not know yet that
    // their global ID may have changed !!!
    // This happens close to intersections of 3+ processes, when a given simplex
    // from process A have one vertex AB shared with B and another AC with C, then C's 
    // ghost simplex that contains AB as a ghost simplex needs update.
    // FIXME: we don t need the fulblown MpiExchg_SimplicesGlobalIds structure
    // here, only the ghost vertices global ids ...
    // FIXME: Is there a better way to do this than yet another MPI call ???
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];
	for (unsigned long j=0;j<sendGlobalIds[i].size();j++)
	  sendGlobalIds[i][j].set(ghostExchange.send[src][j],myRank);
      }

    ghostExchange.exchangeStruct(sendGlobalIds,receivedGlobalIds,
				 mpiType_simplicesGlobalIds.getType(),true,
				 mpiTagsStart_repart+7);

#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int src = ghostExchange.receiveRank[i];
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    GhostSimplex *g=ghostExchange.receive[src][j];	   
	    for (int k=0;k<Simplex::NVERT;k++)
	      {
		Vertex *v=g->getVertex(k);
		if (v->isGhost())
		  {
#pragma omp critical
		    v->setGlobalIdentity(receivedGlobalIds[i][j].vertices[k]);
		  }
	      }	
	  }
      }
    /*
    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
      {
	int haveShadowNei=false;
	int gvCount=0;
	GhostSimplex *s=*it;

	for (int i=0;i<Simplex::NNEI;++i)
	  {	    
	    if (s->getNeighbor(i)->isShadow()) haveShadowNei=true;
	    if (s->getVertex(i)->isGhost()) gvCount++;
	  }

	if ((haveShadowNei)&&(gvCount<NDIM))
	  {
	    glb::console->print<LOG_STD_ALL>
	      ("WARNING: This ghost simplex may not have a shadow neighbor !\n");
	    static_cast<Simplex*>(s)->template print<LOG_STD_ALL>("Ghost:");
	    for (int n=0;n<Simplex::NVERT;++n)
	      s->getVertex(n)->template print<LOG_STD_ALL>();
	  }
      }
    */
    // Finally, we still need to correct the shadows !
    // We start by sending a query for each ghost that has a shadow neighbor (there is 
    // at most one for each ghost). Note that we send the query from the ghost to their
    // orginal copy, so we must reverse the senders/receivers in ghostExchange
    std::vector< std::vector<MpiExchg_ShadowSimplicesFromGhostQuery> > 
      sendSSFGQ(ghostExchange.receiveRank.size());
    std::vector< std::vector<MpiExchg_ShadowSimplicesFromGhostQuery> > 
      receivedSSFGQ;
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int src = ghostExchange.receiveRank[i];
	sendSSFGQ[i].reserve(ghostExchange.receive[src].size()/2);
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    MpiExchg_ShadowSimplicesFromGhostQuery ssfg(ghostExchange.receive[src][j],j);
	    if (!ssfg.isEmpty()) sendSSFGQ[i].push_back(ssfg);
	  }		
      }
    
    ghostExchange.reversedExchangeStruct(sendSSFGQ,receivedSSFGQ,
					 mpiType_shadowSimplicesFromGhostQuery.getType(),
					 true,
					 mpiTagsStart_repart+8);
 
    // now that the query was received, retrieve the shadow's information !
    // in particular, we are interested in the shadow vertex's coord and global ID as
    // well as the new global ID of the shadow
    std::vector< std::vector<MpiExchg_ShadowSimplicesFromGhostReply> > 
      sendSSFGR(ghostExchange.sendRank.size());
    std::vector< std::vector<MpiExchg_ShadowSimplicesFromGhostReply> > 
      receivedSSFGR;

#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];	
	sendSSFGR[i].reserve(receivedSSFGQ[i].size());
	for (unsigned long j=0;j<receivedSSFGQ[i].size();j++)
	  {
	    Simplex *ref = ghostExchange.send[src][receivedSSFGQ[i][j].getBaseCellIndex()];

	    MpiExchg_ShadowSimplicesFromGhostReply 
	      ssfg(ref,receivedSSFGQ[i][j].index,myRank);
	    sendSSFGR[i].push_back(ssfg);
	    //ref->template print<LOG_STD_ALL>("ref =");
	  }
      }
  
    ghostExchange.exchangeStruct(sendSSFGR,receivedSSFGR,
				 mpiType_shadowSimplicesFromGhostReply.getType(),true,
				 mpiTagsStart_repart+9);

    // we need to reset ghostVertexHash because the global IDs have changed !
    // FIXME: do that between the Isend/Ireceive and wait ?
    glb::console->print<LOG_DEBUG>("ghostVertexHash (1): size = %ld  => %.2g \n",
				   ghostVertexHash.size(),
				   double(ghostVertexHash.size())/
				   ghostVertexHashLoadGuess);
    ghostVertexHash.clear_no_resize();
    for (ghostVertexPtr_iterator it=LocalMesh::ghostVertexBegin();it!=itgv_end;++it)
      ghostVertexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
          
    VertexDenseHash shadowVertexHash;    
    set_hash_empty_key(shadowVertexHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(shadowVertexHash,GlobalIdentity::max.get());
    const long shadowVertexHashLoadGuess=1.5*LocalMesh::getNShadowVertices();
    shadowVertexHash.rehash(shadowVertexHashLoadGuess/
			    (shadowVertexHash.max_load_factor()-0.01));
 
     // recycle all the former shadowVertices
    const shadowVertexPtr_iterator itsv_end=LocalMesh::shadowVertexEnd();
    for (shadowVertexPtr_iterator it=LocalMesh::shadowVertexBegin();it!=itsv_end;++it)
      LocalMesh::shadowVertexPool.recycle(*it);
   
    // Now we can rebuild the shadows !!!
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
      {
	int src = ghostExchange.receiveRank[i];

	for (unsigned long j=0;j<receivedSSFGR[i].size();++j)
	  {
	    GhostSimplex *ref = ghostExchange.
	      receive[src][sendSSFGQ[i][j].getBaseCellIndex()];

	    ShadowSimplex *shadow = ref->getNeighbor(sendSSFGQ[i][j].index)->getAsShadow();
	    //if (!shadow->isSet()) // first visit, need to set the vertices !
	      {
	    	for (unsigned int k=0;k<static_cast<unsigned int>(Simplex::NVERT);k++)
	    	  {
	    	    GlobalIdentityValue vgid = receivedSSFGR[i][j].vertexGid[k];
	    	    VertexDenseHash_it it = ghostVertexHash.find(vgid);
		    ShadowVertex *v;
	
		    if (it == ghostVertexHash.end()) // new shadow vertex
	    	      {					
			it = shadowVertexHash.find(vgid);
			if (it==shadowVertexHash.end())
			  {
			    if (k!=receivedSSFGR[i][j].vertexIndex) 
			      {				
			     	char tmp[255];
			     	sprintf(tmp,"err @ k=%ld (should be %ld) from:\n",(long)k,
					receivedSSFGR[i][j].vertexIndex);
			     	static_cast<Simplex*>(ref)->
				  template print<LOG_STD_ALL>(tmp);
				receivedSSFGR[i][j].template print<LOG_STD_ALL>();
				if (shadow!=NULL) shadow->template print<LOG_STD_ALL>();
				
			    
				PRINT_SRC_INFO(LOG_ERROR);
				glb::console->print<LOG_ERROR>
				  ("There are inconsistencies in the global identities, vertex (%d,%d) not found in ghostVertexHash\n",
				   GlobalIdentity(vgid).rank(),GlobalIdentity(vgid).id());

				for (int n=0;n<Simplex::NVERT;++n)
				  ref->getVertex(n)->template print<LOG_STD_ALL>();
				
				exit(-1);
			      }
			    LocalMesh::shadowVertexPool.pop(&v);
			    v->set(receivedSSFGR[i][j].vertexPos,
				   LocalMesh::getNShadowVertices()-1);
			    v->setGlobalIdentity(vgid);
			    v->setGeneration(receivedSSFGR[i][j].generation);
			    v->setData(receivedSSFGR[i][j].vData);
			    shadowVertexHash.insert(std::make_pair(vgid,v));
			  }
			else v=static_cast<ShadowVertex*>(it->second);
		      }
	    	    else v=static_cast<ShadowVertex*>(it->second);
		    // if (shadow->getVertex(k) != NULL) 
		    //   {
		    // 	if (shadow->getVertex(k) != v) exit(0);
		    //   }
		    shadow->setVertex(k,v);
	    	  }
	    	shadow->setGlobalIdentity(receivedSSFGR[i][j].simplexGid);
		shadow->setGeneration(receivedSSFGR[i][j].simplexGeneration);
		shadow->setData(receivedSSFGR[i][j].sData);
	    	shadow->setSetF();
		shadow->setSafeF();
	      }
	    // now still have to add the new neighbor
	    shadow->setNeighbor(receivedSSFGR[i][j].vertexIndex,ref);
	  }
      }

    glb::console->print<LOG_DEBUG>("ghostVertexHash (2): size = %ld  => %.2g \n",
				   ghostVertexHash.size(),
				   double(ghostVertexHash.size())/
				   ghostVertexHashLoadGuess);
    glb::console->print<LOG_DEBUG>("shadowVertexHash: size = %ld  => %.2g \n",
				   shadowVertexHash.size(),
				   double(shadowVertexHash.size())/
				   shadowVertexHashLoadGuess);
    ghostVertexHash.clear();
    shadowVertexHash.clear();

    // now rebuild the shadowExchange structure
    // we will be reusing the ghostOrigin and ghostSendStats structures
    std::vector<int> &shadowSendStats = ghostSendStats;
    std::vector<int> &shadowOrigin = ghostOrigin;
    std::fill(shadowOrigin.begin(),shadowOrigin.end(),0);
    for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();it!=itss_end;++it)
      shadowOrigin[it->getGlobalIdentity().rank()]++;
    
    // rebuild the receive part
    shadowExchange.receiveRank.clear();
    for (unsigned long i=0;i<shadowOrigin.size();++i)
      {
	shadowExchange.receive[i].clear();
	if (shadowOrigin[i]>0)
	  {
	    shadowExchange.receive[i].reserve(shadowOrigin[i]);
	    shadowExchange.receiveRank.push_back(i);
	  }
      }
    
    for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();it!=itss_end;++it)
      shadowExchange.receive[it->getGlobalIdentity().rank()].push_back(*it);
    /*
    // sort them so that they are easy to retrieve
    typename ShadowSimplex::cmpPtrLess shadowSimplexGidCmpLess;
#pragma omp parallel for
    for (int i=0;i<shadowExchange.receiveRank.size();++i)
      {
	std::sort(shadowExchange.receive[i].begin(),
		  shadowExchange.receive[i].end(),
		  shadowSimplexGidCmpLess);
      }
    */
    requestsCount=0;
    mpiCom->Alltoall(shadowOrigin,shadowSendStats);
    for (unsigned long i=0;i<shadowSendStats.size();++i)
      if (shadowSendStats[i]>0) requestsCount++;

    // and send info to rebuild the send part
    std::vector< std::vector<GlobalIdentityValue> > 
      shadowGid(shadowExchange.receiveRank.size());

    requests[5].resize(shadowExchange.receiveRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<shadowExchange.receiveRank.size();++i)
      {
	int src=shadowExchange.receiveRank[i];
	const std::vector<ShadowSimplex *> &curSend = shadowExchange.receive[src];
	shadowGid[i].resize(curSend.size());
	for (unsigned long j=0;j<curSend.size();++j)
	  shadowGid[i][j]=curSend[j]->getGlobalIdentity().get();
	mpiCom->Isend(&shadowGid[i][0],shadowGid[i].size(),
		      shadowExchange.receiveRank[i],&requests[5][i],mpiTagsStart_repart+5);
      }

    // now we can allocate the shadowExchange.send
    shadowExchange.sendRank.clear();    
    for (unsigned long i=0;i<shadowSendStats.size();i++)
      {
	shadowExchange.send[i].clear();
	if (shadowSendStats[i]>0)
	  {
	    shadowExchange.sendRank.push_back(i);
	    shadowExchange.send[i].reserve(shadowSendStats[i]);	   
	  }
      }

    // clean-up some allocated variable we do not need anymore

    glb::console->print<LOG_DEBUG>("simplicesHash: size = %ld  => %.2g \n",
				   simplicesHash.size(),
				   double(simplicesHash.size())/simplicesHashLoadGuess);
    glb::console->print<LOG_DEBUG>("verticesHash: size = %ld  => %.2g \n",
				   verticesHash.size(),
				   double(verticesHash.size())/verticesHashLoadGuess);
    simplicesHash.clear();
    verticesHash.clear();   
    
    // get the simplices in order of their local index
    std::vector<Simplex*> simplicesArr = LocalMesh::getSimplicesArray();
    // receive the global ids of the shadow simplices we have to send
    for (int i=0;i<requestsCount;i++)
      {
	int source=mpiCom->Probe(mpiTagsStart_repart+5); 
	std::vector<GlobalIdentityValue> received(shadowSendStats[source]);
	mpiCom->Recv(&received[0],received.size(),source,mpiTagsStart_repart+5);
	for (unsigned long j=0;j<received.size();++j)
	  {
	    Simplex *s = simplicesArr[GlobalIdentity(received[j]).id()];
	    shadowExchange.send[source].push_back(s);
	  }	
      }

    // wait for all pending Isend requests, as locally allocated buffers
    // may not be freed before they complete
    for (unsigned long i=0;i<requests.size();i++)
      for (unsigned long j=0;j<requests[i].size();j++)
	mpiCom->Wait(&requests[i][j]);

    mpiCom->barrier();
   
    // and update shadowExchange ...
    shadowExchange.updateNCum();

    // We are now DONE !!!!
    // clean-up a little bit before leaving.
    // FIXME: need to clean the tags ?
   
    if (glb::console->willPrint<LOG_DEBUG>())
      {
	shadowExchange.template print<LOG_DEBUG>("SHADOW");
	ghostExchange.template print<LOG_DEBUG>("GHOST");
      }

    /*
    if (glb::debug) 
      {
	const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=it_end;++it)
	  {
	    Simplex *s = *it;
	    Simplex *p = s->getPartner();
	    if (p==NULL) continue;
	    Vertex *v1=s->getVertex(s->getSplitIndex());
	    Vertex *v2=p->getVertex(p->getSplitIndex());
	    if (v1!=v2) 
	      {
		glb::console->print<LOG_STD_ALL>("PARTNER ERROR : \n");
		s->template print<LOG_STD_ALL>("S:");
		v1->template print<LOG_STD_ALL>();
		p->template print<LOG_STD_ALL>("P:");
		v2->template print<LOG_STD_ALL>();
		glb::console->print<LOG_STD_ALL>("----------------\n");
	      }
	  }
	mpiCom->barrier();
	exit(0);
      }
    */
    glb::console->print<LOG_INFO>("done.\n");
    LocalMesh::template printIteratorWastedRatio<LOG_PEDANTIC>();
    glb::console->unIndent(); 

    updateCellsCount();

    if (glb::console->willPrint<LOG_INFO>())
      glb::console->print<LOG_INFO>("Done.\n");
    else
      glb::console->print<LOG_STD>("done.\n");
       
      
    if (glb::debug>1) 
      dumpToNDnetwork("after_repart",
		      IO::NDNET_WithShadows | 
		      IO::NDNET_WithGhosts,
		      true);    
    
    if (glb::debug) checkConsistencyAndReport<LOG_ERROR>("repart");
    /*
    const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=it_end;++it)
      {
        Coord c[NDIM]={0};
        for (int i=0;i<Simplex::NVERT;++i)
          {         
            const Coord *vc=(*it)->getVertex(i)->initCoords.getPointer();
            for (int j=0;j<NDIM;++j)
              c[j] += vc[j];
            c[0]+=(*it)->getVertex(i)->getLocalIndex();

            if ((*it)->getVertex(i)->getLocalIndex()>=LocalMesh::getNVertices())
              {
                glb::console->print<LOG_STD_ALL>
		  ("Loclindex is %d > %d\n",
		   (*it)->getVertex(i)->getLocalIndex(),LocalMesh::getNVertices());
              }
          }
      }
    glb::console->print<LOG_STD_ALL>("COUNT : %d\n",LocalMesh::getNVertices());

    std::vector<long> vec(5000,0);
    long count=0;
    for (auto it=LocalMesh::vertexBegin();it!=LocalMesh::vertexEnd();++it,++count)
      {
	vec[it->getLocalIndex()]=1;
      }
    if (count!=LocalMesh::getNVertices()) 
      glb::console->print<LOG_STD_ALL>("Vertices count differs! (%d!=%d)\n",
				       count,LocalMesh::getNVertices());
    for (int i=0;i<count;++i)
      if (vec[i]!=1)  glb::console->print<LOG_STD_ALL>("No local index %d\n",i);
    */
    return true;
  }

  void checkShadowUniqueness(std::vector<Vertex *> &nonUnique)
  {
    std::set<Vertex *> resultSet;
    std::set<Vertex*,typename Vertex::cmpPosLess> test;
    typedef typename std::set<Vertex*,typename Vertex::cmpPosLess>::iterator Test_it;
    const shadowVertexPtr_iterator itsv_end=LocalMesh::shadowVertexEnd();
    for (shadowVertexPtr_iterator it=LocalMesh::shadowVertexBegin();it!=itsv_end;++it)
      {
	std::pair<Test_it,bool> result = test.insert(*it);
	if (result.second==false)
	  {
	    //PRINT_SRC_INFO(LOG_ERROR);
	    (*result.first)->template print<LOG_STD_ALL>();
	    it->template print<LOG_STD_ALL>();
	    resultSet.insert(*it);
	    resultSet.insert(*result.first);
	    //exit(-1);
	  }
      }
    nonUnique.assign(resultSet.begin(),resultSet.end());
    
  }

  /** \brief compute the load imbalance factor
   *  \return the imbalance factor
  */
  double getLoadImbalanceFactor()
  {    
    return loadImbalanceFactor;

    // double avg=double(getGlobalNCells(NDIM)) / mpiCom->size();   
    // return mpiCom->max(double(LocalMesh::getNCells(NDIM)) / avg);

    /*
    double res = Tree::getLoadImbalanceFactor();
    if (res<0) return loadImbalanceFactor;
    return res;
    */
  }

  /** \brief Adaptively refine the mesh by splitting segments where needed.
   *   
   * This function calls 
   *  \verbatim
      double solver->checkRefine_getValue(Simplex *s) 
      \endverbatim
   * once for each simplex and ghost simplex in the mesh. For each of them, 
   * checkRefine_getValue should return a double precision number representing how much 
   * the simplex needs to be refined (negative or null values stand for no refinement 
   * needed). If the simplex needs refinement, another function 
   * \verbatim
     int checkRefine_getSplitSegmentIndex(Simplex*,double&) 
     \endverbatim
   * is called with parameter a pointer
   * to the simplex to refine and a reference to the score obtained (so that it can be 
   * updated if needed). This funtion must return the index of the segment to split.
   * In case of refinement conflict (i.e. when several edges of a single simplex 
   * have to be split), simplices with higher score are refined
   * first while others are cancelled for the current pass. The function is then called 
   * recursively and as many passes as needed are executed so that all simplices 
   * have a score <= 0 when the function returns.
   *
   * \param solver a pointer to an object implementing a function 
   * checkRefine_getValue and checkRefine_getSplitSegmentIndex. See description hereabove
   * for more information.
   * \param nThreads how many openMP threads to use in parallel.
   * \param nPassMax The maximum number of passes to attempt before giving up. Set nPassMax
   * to 0 for unbounded number of passes.
   */
  template <class S>
  long refine(S *solver, int nThreads=glb::num_omp_threads, int nPassMax=200)
  {
    static int staticCounter=0;
    int pass=1;
    long nRefined=0;
    long nRefinedLocalTotal=0;    
    long nRefinedTotal=0;    
    const int myRank=mpiCom->rank();
    typename TimerPool::Timer timer;
    timer.start();

    std::vector<char> check;    
    typedef std::pair<double,SegmentHandle> CandidateSegment;
    
    if (nThreads<1) nThreads=glb::num_omp_threads;
    std::vector<std::vector<Simplex *> > cachedCandidateSimplices(nThreads);
    std::vector<std::vector<GhostSimplex *> > cachedCandidateGSimplices(nThreads);
    bool useCachedCandidates=false;

    if (glb::debugMeshRefine)
      {
	char fName[256];
	sprintf(fName,"meshRefineCheck_%5.5d_%5.5d",staticCounter++,mpiCom->rank());
	dumpToNDnetwork(fName,
			IO::NDNET_WithShadows|
			IO::NDNET_WithGhosts|
			IO::NDNET_WithNeighbors);	
	mpiCom->barrier();
      }    

    if ((!glb::console->willPrint<LOG_INFO>())&&(glb::console->willPrint<LOG_STD>())) 
      glb::console->printFlush<LOG_STD>("Refining ... ");   

    do {
      glb::console->printFlush<LOG_INFO>("Refining (P%d) ... ",pass++);
      glb::console->printFlush<LOG_PEDANTIC>("(eval) ");
      
      // First check which simplices need refinement and set cache accordingly
      checkRefine_setSimplicesCache(solver,check,nThreads,
				    cachedCandidateSimplices,
				    cachedCandidateGSimplices,
				    useCachedCandidates);

      glb::console->printFlush<LOG_PEDANTIC>("(cflct1) ");
      
      typedef std::vector<CandidateSegment> CandidateVector;
      typedef typename CandidateVector::iterator candidate_iterator;
      //CandidateVector toRefineThread[nThreads];
      //std::vector<Simplex*> cancelledThread[nThreads];
      std::vector<CandidateVector> toRefineThread(nThreads);
      std::vector< std::vector<Simplex*> > cancelledThread(nThreads);
      
      // Build a list of segments to refine.
      // Also ensure that any simplex that needs refinement is not also refined
      // by another simplex
      // In case of conflict, keep the edge with highest score.   
#pragma omp parallel for num_threads(nThreads)
      for (long th=0;th<nThreads;th++)
	{
	  if (false)//useCachedCandidates)
	    {
	      //FIXME: This does not work
	      const auto it_end=cachedCandidateSimplices[th].end();
	      for (auto it=cachedCandidateSimplices[th].begin();it!=it_end;++it)
		checkRefine_checkConflicts(*it,cancelledThread[th],toRefineThread[th]);
	    }
	  else
	    {
	      const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	      for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,nThreads);
		   it!=it_end;++it)
		checkRefine_checkConflicts(*it,cancelledThread[th],toRefineThread[th]);
	    }

	  /*
	  const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	  for (simplexPtr_iterator it=LocalMesh::simplexBegin(i,nThreads);
	       it!=it_end;++it)
	    {
	      float score = it->cache.pfi.f;
	      if (score<=0) continue;
	      int segIndex = it->cache.pfi.i;
	      // printf("score is %e for simplex %ld (seg %d).\n",
	      // 	     score,(long)(*it)->getLocalIndex(),segIndex);
	      	      
	      SegmentHandle h=it->getSegmentHandle(segIndex);
	      segment_circulator ci_end=h->getCirculator();
	      segment_circulator ci=ci_end;
	      
	      // cycle over all simplices containing segment h
	      while((++ci) != ci_end)
		{
		  bool remove=false;
		    
		  //if (fabs(ci->cache.pfi.f-score) <= score*1.E-5)
		  if (ci->cache.pfi.f == score)
		    {	
		      // This happens when two simplices want to refine the same simplex
		      // with the exact same score. We define an arbitrary but deterministic
		      // order for them, based on their global identity
		      // NB: the order must be consistent over different processes
		      
		      SegmentHandle otherSeg = ci->getSegmentHandle(ci->cache.pfi.i);
		      if (LocalMesh::compareSegmentHandlesLess(h,otherSeg,true))
			remove=true;
		      else
			cancelledThread[i].push_back(*ci);
		     
		    }
		  else if (ci->cache.pfi.f > score)
		    {
		      // A simplex with higher score already want to refine an edge
		      // of this simplex
		      remove=true;		      
		    }
		  else
		    {
		      // A simplex with lower score already want to refine an edge
		      // of this simplex
		      cancelledThread[i].push_back(*ci);		      
		    }

		  if (remove)
		    {	
		      score=0;
		      break;
		    }
		};
	      
	      // now if the score is non-null, we can consider this simplex for
	      // refinement.
	      if (score>0)
		toRefineThread[i].push_back(CandidateSegment(score,h));
	      else 
		cancelledThread[i].push_back(*it); 
	    }
	  */
	}
      //if (toRefineThread[0].size()>0) testRefine=true;

      // Merge each thread result and clean the cache of cancelled simplices       
      for (unsigned long j=0;j<cancelledThread[0].size();j++) 
	cancelledThread[0][j]->cache.ptr=NULL;

      CandidateVector &toRefine = toRefineThread[0];    
      // Note this is needed only when using openMP ...
      if (nThreads>1)
	{
	  long toRefineCum[nThreads];
	  toRefineCum[0]=toRefineThread[0].size();
	  for (long i=1;i<nThreads;i++) 
	    toRefineCum[i] = toRefineCum[i-1]+toRefineThread[i].size();
	  toRefine.resize(toRefineCum[nThreads-1]);
#pragma omp parallel for num_threads(nThreads)
	  for (long i=0;i<nThreads;i++)
	    {
	      if (i==0)
		{
		  for (long j=1;j<nThreads;j++)
		    for (unsigned long k=0;k<cancelledThread[j].size();k++) 
		      cancelledThread[j][k]->cache.ptr=NULL;
		}
	      else
		{
		  std::copy(toRefineThread[i].begin(),toRefineThread[i].end(),
			    toRefine.begin()+toRefineCum[i-1]);
		}
	    }
	}
      
      for (candidate_iterator it=toRefine.begin();it!=toRefine.end();++it)
	it->first=it->second->getSimplex()->cache.pfi.f;	 
      
      // Now we want to import all the candidate segments that belong to ghosts simplices
      // Note that we already computed those whose vertices are shared vertices, but not
      // the ones with at least one ghost vertex !
      if (LocalMesh::getNGhostSimplices()>0) 
	{	  
	  glb::console->printFlush<LOG_PEDANTIC>("(ghosts) ");	 
	  
	  std::vector< std::vector<MpiExchg_RefineQueryResult> > 
	    toSend(ghostExchange.sendRank.size());

	  // This is for testing purpose
	  //std::vector< std::vector< typename MpiExchg_RefineQueryResult::MpiStruct > > 
	  // tmpS(ghostExchange.sendRank.size()); 
	  //std::vector< std::vector< typename MpiExchg_RefineQueryResult::MpiStruct > > 
	  // tmpR; 
	  //std::vector<MPI_Request> reqs; 
	  //ghostExchange.iReceive(tmpR,reqs,mpiType_refineQueryResult); 

#pragma omp parallel for num_threads(nThreads)
	  for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
	    {
	      std::vector<Simplex *> &curSend= 
		ghostExchange.send[ghostExchange.sendRank[i]];

	      for (unsigned long j=0;j<curSend.size();j++)
		{
		  Simplex *s=curSend[j];		
		  
		  if (s->cache.pfi.f>0)	
		    {
		      MpiExchg_RefineQueryResult res;
		      //SegmentHandle h=s->getSegmentHandle(s->cache.pfi.i);
		      res.set(s->cache.pfi.f,s->getSegmentHandle(s->cache.pfi.i));
		      res.setBaseCellIndex(j);
		      toSend[i].push_back(res);			  
		    }		
		}
	      //ghostExchange.iSend(i,toSend[i],tmpS[i],reqs,mpiType_refineQueryResult); //TEST
	    }

	  //std::vector< MpiExchg_RefineQueryResult > received;
	  //ghostExchange.wait(tmpR,received,reqs,mpiType_refineQueryResult);
	
	  std::vector< MpiExchg_RefineQueryResult > received;
	  ghostExchange.exchangeStruct(toSend,received,mpiType_refineQueryResult.getType(),
				       mpiTagsStart_refine+0);	  

	  toRefine.reserve(toRefine.size() + received.size());
	  for (unsigned long i=0;i<received.size();i++)
	    toRefine.push_back(CandidateSegment(received[i].value,received[i].handle));
	}     
      
      glb::console->printFlush<LOG_PEDANTIC>("(cflct2) ");
      
      // CHECKME: this should be faster than resetSimplicesCache but relies
      // on the fact that non split simplices have a cache.ptr==NULL !!!
      for (candidate_iterator it=toRefine.begin();it!=toRefine.end();++it)
	it->second->getSimplex()->cache.ptr=NULL;
      
      // Check !
      // const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
      // for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=it_end;++it)
      // 	if (it->cache.ptr!=NULL) exit(0);
	   
      //LocalMesh::resetSimplicesCache(nThreads);  
      LocalMesh::resetGhostSimplicesCache(nThreads);     
      std::vector<Simplex*> conflictSimplices;
      
      // Now we can check for second order conflicts (i.e a simplex that is
      // refined by at least two other simplices)
      // => for any given simplex, only one edge may be refined at a time
      // After this loop, any cancelled candidate score will be negated
      // FIXME: any safe way to use openMP here ?
      for (candidate_iterator it=toRefine.begin();it!=toRefine.end();++it)
	{
	  // Those ones where cancelled by direct neighbors 
	  // We are being conservative here : those with negative score, that were already 
	  // cancelled by remote conflicts, can still cancel other refinements so that the 
	  // result is deterministic !
	  if (it->first==0) continue;
	  	  
	  segment_circulator ci_end=it->second->getCirculator();
	  segment_circulator ci=ci_end;	  
	  
	  // for each candidate, tag any simplex its refinement will affect
	  // with a pointer to it (in ci->cache.ptr)
	  do
	    {
	      if (ci->cache.ptr==NULL)
		{
		  // this simplex is not yet affected by anything
		  ci->cache.ptr = static_cast<void*>(&(*it));
		}
	      else
		{
		  // this simplex is already affected by a refinement !
		  CandidateSegment *other = static_cast<CandidateSegment*>(ci->cache.ptr);
		  Simplex *otherS=other->second->getSimplex();
		  Simplex *curS=it->second->getSimplex();
		  float diff = fabs(other->first) - fabs(it->first);
		  if (diff>0)
		    {
		      // cancel current candidate
		      curS->cache.ptr = static_cast<void*>(&(*other));
		      it->first = -fabs(it->first);
		    }
		  else if (diff<0)
		    {
		      // cancel the other simplex's candidate
		      ci->cache.ptr=static_cast<void*>(&(*it));	
		      other->first = -fabs(other->first);
		      //otherS->cache.ptr=static_cast<void*>(&(*it));
		    }
		  else
		    {		     
		      // Two candidates with the exact same score are in conflict,
		      GlobalIdentity otherGid = otherS->getGlobalIdentity(myRank);
		      GlobalIdentity curGid = curS->getGlobalIdentity(myRank);
		      if (curGid != otherGid)
			{
			  bool cmp = LocalMesh::compareSegmentHandlesLess
			    (it->second,other->second,true);
			  /*
			  // importance is decided by comparing the simplices global ID
			  // This works but does not give nice results !
			  bool cmp = (otherGid>curGid);
			  */
			  if (cmp)
			    {
			      // cancel current candidate
			      curS->cache.ptr = static_cast<void*>(&(*other));
			      it->first = -fabs(it->first);
			    }
			  else 
			    {
			      // cancel the other simplex's candidate
			      ci->cache.ptr=static_cast<void*>(&(*it));
			      other->first = -fabs(other->first);
			      //otherS->cache.ptr=static_cast<void*>(&(*it));
			    }						  
			}

		      
		      // NB : we don't do anything when the two candidates are the same !
		      // This WILL happen, as we also test previously cancelled candidates
		      // for possible conflicts !
		    }		 
		}
	    } while ((++ci)!=ci_end);
	}

      // Now it is still possible that a ghost simplex was cancelled remotely but
      // not locally. 
      long nRemoteConflicts=0;        
      // Check possible leftover conflicts with remote processes
      if (LocalMesh::getNGhostSimplices()>0) 
	{
	  // Set the cache of every simplex candidate to a pointer to its segment and
	  // score. Note that we DO include the previously cancelled ones (negative score)
	  // but not the initially cancelled ones (score == 0)
	  for (candidate_iterator it=toRefine.begin();it!=toRefine.end();++it)
	    if (it->first != 0) 
	      it->second->getSimplex()->cache.ptr=&(*it);
	  
	  glb::console->printFlush<LOG_PEDANTIC>("(ghosts) ");	 
	  
	  std::vector< std::vector< MpiExchg_RefineQueryResult > > 
	    toSend(ghostExchange.sendRank.size());
	  std::vector< std::vector< MpiExchg_RefineQueryResult > > 
	    toSendR(ghostExchange.receiveRank.size());

	  // Start by checking that a locally non cancelled simplex was not cancelled
	  // on a remote process: check the cancelled ghost simplices, and send the info
	  // to their local process that they should be cancelled. Cancelled candidates 
	  // can be recognized by the fact that their simplices S cache.ptr points to 
	  // a candidate whose simplex is S itself, moreover, the candidate score 
	  // is negative as it was cancelled
#pragma omp parallel for num_threads(nThreads)
	  for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
	    {
	      std::vector<GhostSimplex *> &curSend = 
		ghostExchange.receive[ghostExchange.receiveRank[i]];

	      for (unsigned long j=0;j<curSend.size();j++)
		{
		  Simplex *s=curSend[j];				  
		  if (s->cache.ptr!=NULL) 
		    {
		      CandidateSegment *r = static_cast<CandidateSegment*>(s->cache.ptr);

		      if ((r->second->getSimplex() == s)&&(r->first<0))
			{
			  // This segment is a cancelled candidate
			  MpiExchg_RefineQueryResult res;
			  res.set(r->first,r->second);
			  res.setBaseCellIndex(j);
			  toSendR[i].push_back(res);			 
			}
		    }		
		}
	    }	 

	  std::vector< MpiExchg_RefineQueryResult > receivedR;
	  ghostExchange.reversedExchangeStruct(toSendR,receivedR,
					       mpiType_refineQueryResult.getType(),
					       mpiTagsStart_refine+2);	

	  // compare remote and local cancellations
 	  for (unsigned long i=0;i<receivedR.size();i++)
	    {
	      Simplex *s=receivedR[i].handle->getSimplex();	     
	      // any simplex considered here was cancelled remotely
	      if (s->cache.ptr!=NULL) 
		{		  
		  // but this one was NOT cancelled locally !
		  nRemoteConflicts++;	
		  //Remote cancellation -> set negative value to cancel refinement
		  CandidateSegment *r = static_cast<CandidateSegment*>(s->cache.ptr);
		  r->first=-1;
		}
	    }

	  // Now communicate the locally cancelled simplices to other processes ghosts so
	  // that everything is consistent
#pragma omp parallel for num_threads(nThreads)
	  for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
	    {
	      std::vector<Simplex *> &curSend = 
		ghostExchange.send[ghostExchange.sendRank[i]];

	      for (unsigned long j=0;j<curSend.size();j++)
		{
		  Simplex *s=curSend[j];				  
		  if (s->cache.ptr!=NULL) 
		    {
		      CandidateSegment *r = static_cast<CandidateSegment*>(s->cache.ptr);
		      if ((r->second->getSimplex() == s)&&(r->first<0))
			{			  
			  // This segment is a cancelled candidate
			  MpiExchg_RefineQueryResult res;
			  res.set(r->first,r->second);
			  res.setBaseCellIndex(j);
			  toSend[i].push_back(res);			 
			}
		    }		
		}
	    }

	  std::vector< MpiExchg_RefineQueryResult > received;
	  ghostExchange.exchangeStruct(toSend,received,
				       mpiType_refineQueryResult.getType(),
				       mpiTagsStart_refine+1);	

	  // compare remote and local cancellations
 	  for (unsigned long i=0;i<received.size();i++)
	    {
	      Simplex *s=received[i].handle->getSimplex();	     
	      // any simplex considered here was cancelled remotely
	      if (s->cache.ptr!=NULL) 
		{
		  // but this one was NOT cancelled locally !
		  nRemoteConflicts++;	
		  // set cache to NULL so that we know it should be cancelled
		  s->cache.ptr=NULL; 
		}
	    }
	}

      // Enable cached candidates for next pass
      useCachedCandidates=true;

      long nGhostRefined=0;
      unsigned long toRefineCount=0;
      // We can now build a final list of segments to refine ...
      for (unsigned long i=0;i<toRefine.size();++i)
	{
	  Simplex *s=toRefine[i].second->getSimplex();
	  if (toRefine[i].first<=0) 	  
	    {
	      s->cache.ptr=NULL;
	      continue;
	    }	  
	  // this happens when a segment refinement is remotely cancelled
	  if (s->cache.ptr==NULL) continue;		  
	  if (s->isShadowOrGhost()) nGhostRefined++;
	  if (i!=toRefineCount) std::swap(toRefine[toRefineCount],toRefine[i]);
	  //toRefine[toRefineCount]=toRefine[i];
	  toRefineCount++;		  
	}

      // Clean up cache before next pass
      /*
      if (useCachedCandidates)
	{
	  unsigned long sz=toRefine.size();
#pragma omp parallel for num_threads(nThreads)
	  for (long i=0;i<sz;++i)
	    {
	      segment_circulator ci_end=toRefine[i].second->getCirculator();
	      segment_circulator ci=ci_end;
	      do
		{
		  ci->cache.ptr=NULL;
		} while ((++ci)!=ci_end);
	    }
	}
      */

      // Keep only the simplices we will really refine !
      toRefine.resize(toRefineCount);
      nRefined=toRefine.size()-nGhostRefined;

      /*
       for (candidate_iterator it=toRefine.begin();it!=toRefine.end();++it)
	{
	  if ((((*it).second->getVertex(0)->getGlobalIdentity()==GlobalIdentity(0,258))||
	       ((*it).second->getVertex(1)->getGlobalIdentity()==GlobalIdentity(0,258)))&&
	      (((*it).second->getVertex(0)->getGlobalIdentity()==GlobalIdentity(0,381))||
	       ((*it).second->getVertex(1)->getGlobalIdentity()==GlobalIdentity(0,381))))
	    {
	      (*it).second->template print<LOG_STD_ALL>();
	      (*it).second->getVertex(0)->template print<LOG_STD_ALL>();
	      (*it).second->getVertex(1)->template print<LOG_STD_ALL>();
	      (*it).second->getSimplex()->template print<LOG_STD_ALL>("AFSIMPLEX:");
	      
	    }
	}
      */

      // Let's split all the local simplices that need to be 
      glb::console->printFlush<LOG_PEDANTIC>("(refine) "); 
      
      UMapGlobal newSharedVerticesMap;
      std::vector<Vertex*> newVertices;
      std::vector<Simplex*> newSimplices;
      std::vector<unsigned long> nNewSimplicesCum;  
      /*
      if (toRefine.size()<=6)
	{
	  glb::console->printFlush<LOG_STD>("\n");
	  for (int i=0;i<toRefine.size();++i)
	    {
	      Simplex *s=toRefine[i].second->getSimplex();
	      glb::console->print<LOG_STD>("Refining simplex %ld:\n",
					   (long)s->getLocalIndex());
	      
	      
	      s->template print<LOG_STD>("REF");
	      for (int j=0;j<Simplex::NVERT;++j)
		{
		  Vertex *v=s->getVertex(j);
		  v->template print<LOG_STD>();
		}
	      std::cout << "(" << s->segTracers.getValueAt(0);
	      for (int j=1;j<Simplex::NSEG*NDIM_W;++j) 
		{
		  if ((j%NDIM_W) == 0) 
		    std::cout << "(" << s->segTracers.getValueAt(j);
		  else
		    std::cout << "," << s->segTracers.getValueAt(j);
		  if ((j%NDIM_W) == NDIM_W-1) std::cout << ")\n";
		}
	
	    }
	}
      */
      // NOTE : after calling splitSegments, a simplex S will be split into 
      // two simplices A and B. Only one new simplex is actually created, B, and
      // A is only a shrunk S (i.e. A will be stored in memory where S was).
      // For convenience, the cache.ptr of A will point to the simplex that made S 
      // split, while the cache.ptr of B will point to the newly created vertex.
      // It is also easy to differentiate A from B as A will get a TAG flag 
      // (setTaggedF()) while B won't. 
      LocalMesh::splitSegments(toRefine,newSharedVerticesMap,
      			       newVertices,newSimplices,
			       nNewSimplicesCum,
			       solver,nThreads);
      
      // Add split simplices to  cachedCandidateSimplices for next pass
#pragma omp parallel for num_threads(nThreads)
      for (int th=0; th<nThreads; th++)
	{
	  long delta=newSimplices.size()/nThreads;
	  long start = th*delta;
	  long stop = (th+1)*delta;
	  if (th==nThreads-1) stop=newSimplices.size();
	  for (long j=start;j<stop;++j)
	    {
	      cachedCandidateSimplices[th].push_back(newSimplices[j]);
	      cachedCandidateSimplices[th].push_back(newSimplices[j]->getPartner());
	    }
	}      
      
      // now that we are split, we still have to synchronize the ghosts/shadows
      // properties (e.g. global IDs, neighbors, ghostExchange structures, ...)
      glb::console->printFlush<LOG_PEDANTIC>("(sync) ");	
             
      std::vector<Vertex*> newVerticesG;
      std::vector<Simplex*> newSimplicesG;
      std::vector<unsigned long> nNewSimplicesCumG;

      if (LocalMesh::getNGhostSimplices()>0)
	{
	  synchronizeAfterSplitting(newSharedVerticesMap,
				    LocalMesh::ghostVertexPool, 
				    LocalMesh::ghostSimplexPool,
				    ghostExchange,newVerticesG,
				    newSimplicesG,nNewSimplicesCumG,nThreads);
	  
	   // Add split ghost simplices to  cachedCandidateGSimplices for next pass
#pragma omp parallel for num_threads(nThreads)
	  for (int th=0; th<nThreads; th++)
	    {
	      long delta=newSimplicesG.size()/nThreads;
	      long start = th*delta;
	      long stop = (th+1)*delta;
	      if (th==nThreads-1) stop=newSimplicesG.size();
	      for (long j=start;j<stop;++j)
		{
		  // cache.c[0] is the index of the neighbor that contains the 
		  // split partner
		  GhostSimplex *newSimplex=
		    newSimplicesG[j]->getAsGhost();
		  GhostSimplex *otherSimplex=
		    newSimplex->getNeighbor(newSimplex->cache.c[0])->getAsGhost();
		  
		  cachedCandidateGSimplices[th].push_back(newSimplex);
		  cachedCandidateGSimplices[th].push_back(otherSimplex);
		}
	    }
	  
	}

      std::vector<Vertex*> newVerticesS;
      std::vector<Simplex*> newSimplicesS;
      std::vector<unsigned long> nNewSimplicesCumS;
      
      if (LocalMesh::getNShadowSimplices()>0)
	synchronizeAfterSplitting(newSharedVerticesMap,
				  LocalMesh::shadowVertexPool, 
				  LocalMesh::shadowSimplexPool,
				  shadowExchange,newVerticesS,
				  newSimplicesS,nNewSimplicesCumS,nThreads);
      
      // We can now rebuild the neighborhood of all the split simplices
      // All newSimplices/newVertices must be created before calling this
      LocalMesh::fixSimplicesNeighborsAfterSplitting(newVertices,newSimplices,
						     nNewSimplicesCum,false);
      if (LocalMesh::getNGhostSimplices()>0)
	LocalMesh::fixSimplicesNeighborsAfterSplitting(newVerticesG,newSimplicesG,
						       nNewSimplicesCumG,false);
      if (LocalMesh::getNShadowSimplices()>0)
	LocalMesh::fixSimplicesNeighborsAfterSplitting(newVerticesS,newSimplicesS,
						       nNewSimplicesCumS,false);  
      /*
      // Need to clean the cache if we are going to use cached candidates
      if (useCachedCandidates)
	{
	  unsigned long sz=toRefine.size();
#pragma omp parallel for num_threads(nThreads)
	  for (long i=0;i<sz;++i)
	    {
	      if (toRefine[i].first==0) continue;
	      segment_circulator ci_end=toRefine[i].second->getCirculator();
	      segment_circulator ci=ci_end;
	      do
		{
		  ci->cache.ptr=NULL;
		} while ((++ci)!=ci_end);
	    }
	}      
      */

      /*
      // Add split simplices to  cachedCandidateSimplices for next pass
#pragma omp parallel for num_threads(nThreads)
      for (int th=0; th<nThreads; th++)
	{
	  long delta=newSimplices.size()/nThreads;
	  long start = th*delta;
	  long stop = (th+1)*delta;
	  if (th==nThreads-1) stop=newSimplices.size();
	  for (long j=start;j<stop;++j)
	    {
	      // The two partner split simplices
	      Simplex *s1=newSimplices[j];
	      Simplex *s2=newSimplices[j]->getNeighbor(s1->cache.c[0]);	   
	      Vertex *vSplit=s1->getVertex(s1->cache.c[1]);
	      // Any simplex incident to an edge of the new simplex not containing
	      // the split vertex
	      for (int i=0;i<Simplex::NSEG;++i)
		{
		  auto sh=s1->getSegmentHandle(i);
		  if ((sh->getVertex(0)!=vSplit)&&(sh->getVertex(1)!=vSplit))
		    sh->getAdjacentSimplices
		      (std::back_inserter(cachedCandidateSimplices[th]));
		}

	      // Any simplex incident to an edge of the new simplex partner not containing
	      // the split vertex and that is not shared with the new simplex
	      Vertex *oppVertex=s2->getVertex(s2->cache.c[0]);
	      for (int i=0;i<Simplex::NSEG;++i)
		{
		  auto sh=s2->getSegmentHandle(i);
		  Vertex *v0=sh->getVertex(0);
		  Vertex *v1=sh->getVertex(1);
		  if ((v0!=vSplit)&&
		      (v1!=vSplit)&&
		      ((v0==oppVertex)||(v1==oppVertex)))
		    sh->getAdjacentSimplices
		      (std::back_inserter(cachedCandidateSimplices[th]));
		}	     
	    }	 
	}
      */
      long nRefinedLocal=nRefined;
      nRefinedLocalTotal+=nRefined; 
     
      // FIXME: Could we avoid the GLOBAL synchronization here ????
      // probably not as we need to know whether another pass is needed ...
      nRefined=mpiCom->sum(nRefined); 
      nRefinedTotal+=nRefined;
      
      if (glb::debug) checkConsistencyAndReport<LOG_ERROR>("refine");

      if (nRemoteConflicts>0)
	{
	  if (nRefined)
	    glb::console->print<LOG_INFO>
	      ("done. (%ld segments refined, %ld remote conflicts)\n",
	       nRefined,nRemoteConflicts);
	  else
	    glb::console->print<LOG_INFO>("done. (%ld remote conflicts)\n",
					  nRemoteConflicts);
	}
      else
	{
	  if (nRefined)
	    glb::console->print<LOG_INFO>("done. (%ld/%ld segments refined)\n",
					  nRefinedLocal,nRefined);
	  else
	    glb::console->print<LOG_INFO>("done.\n");
	}
      
      if ((nPassMax>0)&&(pass>nPassMax)&&(nRefined>0))
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Too many iterations, giving up.");
	  glb::console->print<LOG_ERROR>("Check your refinement criterion ...\n");
	  exit(-1);
	}

    } while ((nRefined>0)&&(!glb::debugMeshRefine));
 
    updateCellsCount();
    double timeSpent=timer.stop();    
    
    if ((!glb::console->willPrint<LOG_INFO>())&&(glb::console->willPrint<LOG_STD>())) 
      glb::console->print<LOG_STD>("done in %.2fs.\n",timeSpent);   
    else if (nRefinedTotal>0)
      glb::console->print<LOG_INFO>
	("A total of %ld/%ld segments were Locally/Globally refined in %.2fs.\n",
	 nRefinedLocalTotal,nRefinedTotal,timeSpent);
    else
      glb::console->print<LOG_INFO>
	("No refinement needed (checking took %.2fs).\n",timeSpent);
    
    //mpiCom->barrier();
    
    //if (testRefine) glb::console->print<LOG_STD_ALL>(" THERE ARE LEFTOVERS !!!!\n");  
    //mpiCom->barrier();
    //LocalMesh::template printIteratorWastedRatio<LOG_PEDANTIC_ALL>();

    if (glb::debugMeshRefine)
      {
	//static int index=1;
	if (nRefined==0)
	  {
	    char fName[256];
	    sprintf(fName,"meshRefineCheck_%5.5d_%5.5d",staticCounter++,mpiCom->rank());
	    dumpToNDnetwork(fName,
			    IO::NDNET_WithShadows| 
			    IO::NDNET_WithGhosts|
			    IO::NDNET_WithNeighbors,
			    true);
	    mpiCom->barrier();
	  }
	else nRefinedLocalTotal += refine(solver,nThreads);
      }

    return nRefinedLocalTotal;
  }

  /** 
   * \brief Simulates a simplex splitting along every possible edge, without 
   * affecting the actual mesh. Note that vertices' boundary flag and simplices non-split 
   * neighbors are not updated.
   * \param s    the simplex to split
   * \param[out] newVertex An array of NSEG new vertices, one for each split segment
   * \param[out] simplices A vector of the simplices incident to the split segments
   * \param[out] splitSimplices An array of vectors containing the split simplices. 
   *  splitSimplices[seg][0/1][n] will contain the two simplices generated by breaking 
   * segment \a seg of the nth simplex adjacent to segment \a seg of simplex \a s. 
   * More specifically, \a splitSimplices[seg][0][n] contains the simplex generated by 
   * shrinking simplex \a simplices[seg][n] while \a splitSimplices[seg][1][n] contains 
   * its newly created partner. If 'splitNeighbors' is false, then simplices[seg] has 
   * only one element (simplices[seg][0]=s) and it is split into the two simplices 
   * splitSimplices[seg][0/1][0].
   * \param updateVertexData Whether to update the vertex data (may save time if updated
   * data is not needed).
   * \param updateSimplexData Whether to update the simplex data (may save time if updated
   * data is not needed).
   * \param splitNeighbors if false, only simplex 's' is actually split, otherwhise all 
   * simplices adjacent to segments of \a s are split.
   * \tparam UpdateDataIndex If UpdateDataIndex>=0 and updateSimplexData=true, then 
   * update only the simplex data field with index 'UpdateDataIndex'. Can be used
   * to save to time ...
   * \warning Vertices' boundary flag and simplices non-split 
   * neighbors are not updated!
   */
  template <int UpdateDataIndex=-1>
    void simulateAllSimplexSplits(Simplex *s, 
				  Vertex newVertex[Simplex::NSEG],
				  std::vector<Simplex*> 
				  simplices[Simplex::NSEG],
				  std::vector<Simplex> 
				  splitSimplices[Simplex::NSEG][2],
				  bool updateVertexData=true,
				  bool updateSimplexData=true,
				  bool splitNeighbors=true) const
  {             
    for (int seg=0;seg<Simplex::NSEG;++seg)
      {
	simulateSimplexSplit<UpdateDataIndex>(s,seg,
					      newVertex[seg],
					      simplices[seg],
					      &splitSimplices[seg][0],
					      updateVertexData,
					      updateSimplexData,
					      splitNeighbors);
      }
    /*
    //std::vector<Simplex*> sList;
    const long rBufferSize=Vertex::template 
      getSimplexRefineBufferSize<MyType,SegmentHandle>();

    for (int seg=0;seg<Simplex::NSEG;++seg)
      {
	std::vector<Simplex*> &sList = simplices[seg]; 
	SegmentHandle handle=s->getSegmentHandle(seg);
	sList.clear();
	if (splitNeighbors)
	  handle->getAdjacentSimplices(std::back_inserter(sList));
	else
	  sList.push_back(s);
	char buffer[rBufferSize*sList.size()];

	// Create the new vertex
	if (updateVertexData)
	  {
	    // Update coordinates and data
	    Vertex::Data::onRefineVertex
	      (this,handle,&newVertex[seg],&sList[0],sList.size(),buffer);
	  }
	else
	  {
	    // Only the coordinates of the vertex are updated in that case
	    newVertex[seg].setRefinedCoords(this,handle,&sList[0],sList.size(),buffer);
	  }	  
	  
	if (splitNeighbors) 
	  {	  
	    splitSimplices[seg][0].resize(sList.size());
	    splitSimplices[seg][1].resize(sList.size());
	    for (long n=0;n<sList.size();++n)
	      splitSimplices[seg][0][n]=(*sList[n]);	      
	  }
	else
	  {
	    splitSimplices[seg][0].resize(1);
	    splitSimplices[seg][1].resize(1);
	    splitSimplices[seg][0][0]=(*s);
	  }

	char *curBuffer = buffer;
	// split simplices adjacent to the broken segment
	for (int n=0;n<splitSimplices[seg][0].size();++n)
	  {	      
	    splitSimplices[seg][0][n].splitSegment
	      (handle,NULL,&newVertex[seg],&splitSimplices[seg][1][n]);

	    // update simplex data if required
	    if (updateSimplexData)
	      {
		Simplex *s0=&splitSimplices[seg][0][n];
		Simplex *s1=&splitSimplices[seg][1][n];	   	      
	      
		if (UpdateDataIndex>=0)
		  {
		    // Update only data with index UpdateDataIndex
		    Simplex::Data::template refineSimplexData<UpdateDataIndex>
		      (this,handle,&newVertex[seg],&s0,&s1,1,curBuffer);
		  }
		else
		  {
		    // update all data
		    Simplex::Data::onRefineSimplices
		      (this,handle,&newVertex[seg],&s0,&s1,1,curBuffer);
		  }
	      }
	    curBuffer += rBufferSize;
	  }
	
	// and fix the split simplices neighbors
	// for (int n=0;n<splitSimplices[seg][0].size();++n)
	// {	       
	// splitSimplices[seg][0][n].fixNeighborsAfterSplitting(false,false);
	// splitSimplices[seg][1][n].fixNeighborsAfterSplitting(false,false);
	// }		

      }
*/
  }


  /** 
   * \brief Simulates a simplex splitting along edge \a seg, without 
   * affecting the actual mesh. Note that vertices' boundary flag and simplices non-split 
   * neighbors are not updated.
   * \param s    the simplex to split
   * \param seg  the index of the edge of \a s to split
   * \param[out] newVertex The new vertex introduced when spliting edge \a seg
   * \param[out] simplices A vector of the simplices incident to the split segment
   * \param[out] splitSimplices An array of vectors containing the split simplices. 
   *  splitSimplices[0/1][n] will contain the two simplices generated by breaking 
   * segment \a seg of the nth simplex adjacent to segment \a seg of simplex \a s. 
   * More specifically, \a splitSimplices[0][n] contains the simplex generated by 
   * shrinking simplex \a simplices[n] while \a splitSimplices[1][n] contains 
   * its newly created partner. If 'splitNeighbors' is false, then simplices has 
   * only one element (simplices[0]=s) and it is split into the two simplices 
   * splitSimplices[0/1][0].
   * \param updateVertexData Whether to update the vertex data (may save time if updated
   * data is not needed).
   * \param updateSimplexData Whether to update the simplex data (may save time if updated
   * data is not needed).
   * \param splitNeighbors if false, only simplex 's' is actually split, otherwhise all 
   * simplices adjacent to segment \a seg of \a s are split.
   * \tparam UpdateDataIndex If UpdateDataIndex>=0 and updateSimplexData=true, then 
   * update only the simplex data field with index 'UpdateDataIndex'. Can be used
   * to save to time ...
   * \warning Vertices' boundary flag and simplices non-split 
   * neighbors are not updated!
   */
  template <int UpdateDataIndex=-1>
    void simulateSimplexSplit(Simplex *s, int seg,
			      Vertex &newVertex,
			      std::vector<Simplex*> &simplices,
			      std::vector<Simplex> splitSimplices[2],
			      bool updateVertexData=true,
			      bool updateSimplexData=true,
			      bool splitNeighbors=true) const
  {   
     //std::vector<Simplex*> sList;
    const long rBufferSize=Vertex::template 
      getSimplexRefineBufferSize<MyType,SegmentHandle>();

    std::vector<Simplex*> &sList = simplices; 
    SegmentHandle handle=s->getSegmentHandle(seg);
    sList.clear();
    if (splitNeighbors)
      handle->getAdjacentSimplices(std::back_inserter(sList));
    else
      sList.push_back(s);
    char buffer[rBufferSize*sList.size()];

    // Create the new vertex
    if (updateVertexData)
      {
	// Update coordinates and data
	Vertex::Data::onRefineVertex
	  (this,handle,&newVertex,&sList[0],sList.size(),buffer);
      }
    else
      {
	// Only the coordinates of the vertex are updated in that case
	newVertex.setRefinedCoords(this,handle,&sList[0],sList.size(),buffer);
      }	  
	  
    if (splitNeighbors) 
      {	  
	splitSimplices[0].resize(sList.size());
	splitSimplices[1].resize(sList.size());
	for (long n=0;n<sList.size();++n)
	  splitSimplices[0][n]=(*sList[n]);	      
      }
    else
      {
	splitSimplices[0].resize(1);
	splitSimplices[1].resize(1);
	splitSimplices[0][0]=(*s);
      }

    char *curBuffer = buffer;
    // split simplices adjacent to the broken segment
    for (int n=0;n<splitSimplices[0].size();++n)
      {	      
	splitSimplices[0][n].splitSegment
	  (handle,NULL,&newVertex,&splitSimplices[1][n]);

	// update simplex data if required
	if (updateSimplexData)
	  {
	    Simplex *s0=&splitSimplices[0][n];
	    Simplex *s1=&splitSimplices[1][n];	   	      
	      
	    if (UpdateDataIndex>=0)
	      {
		// Update only data with index UpdateDataIndex
		Simplex::Data::template refineSimplexData<UpdateDataIndex>
		  (this,handle,&newVertex,&s0,&s1,1,curBuffer);
	      }
	    else
	      {
		// update all data
		Simplex::Data::onRefineSimplices
		  (this,handle,&newVertex,&s0,&s1,1,curBuffer);
	      }
	  }
	curBuffer += rBufferSize;
      }

    /*
    // and fix the split simplices neighbors
    for (int n=0;n<splitSimplices[seg][0].size();++n)
    {	       
    splitSimplices[seg][0][n].fixNeighborsAfterSplitting(false,false);
    splitSimplices[seg][1][n].fixNeighborsAfterSplitting(false,false);
    }	 
    */
  }


  /** 
   * \brief Simulates a simplex splitting along every possible edge, without 
   * affecting the actual mesh. Only simplex \a s is split, so this is equivalent 
   * to calling simulateAllSimplexSplits with \a splitNeighbors=false, but faster. 
   * Note that vertices' boundary flag and 
   * simplices non-split neighbors are not updated.
   * \param s    the simplex to split
   * \param[out] newVertex An array of NSEG new vertices, one for each split segment
   * \param[out] splitSimplices An array of vectors containing the split simplices. 
   * splitSimplices[seg][0/1] will contain the two simplices generated from \a s by 
   * breaking  edge \a seg. More specifically, \a splitSimplices[seg][0] contains the 
   * simplex generated by shrinking \a s while \a splitSimplices[seg][1] contains 
   * its newly created partner. 
   * \param updateVertexData Whether to update the vertex data (may save time if updated
   * data is not needed).
   * \param updateSimplexData Whether to update the simplex data (may save time if updated
   * data is not needed).   
   * \tparam UpdateDataIndex If UpdateDataIndex>=0 and updateSimplexData=true, then 
   * update only the simplex data field with index 'UpdateDataIndex'. Can be used
   * to save to time ...
   * \warning Vertices' boundary flag and simplices non-split 
   * neighbors are not updated!
   */
  template <int UpdateDataIndex=-1>
    void simulateAllIsolatedSimplexSplits(Simplex *s, 
					  Vertex newVertex[Simplex::NSEG],
					  Simplex splitSimplices[Simplex::NSEG][2],
					  bool updateVertexData=true,
					  bool updateSimplexData=true) const
  {             
    for (int seg=0;seg<Simplex::NSEG;++seg)
      {
	simulateIsolatedSimplexSplit<UpdateDataIndex>(s,seg,
						      newVertex[seg],
						      &splitSimplices[seg][0],
						      updateVertexData,
						      updateSimplexData);
      }
  }

  /** 
   * \brief Simulates a simplex splitting along edge \a seg, without 
   * affecting the actual mesh. Only simplex \a s is split, so this is equivalent 
   * to calling simulateSimplexSplit with \a splitNeighbors=false, but faster. 
   * Note that vertices' boundary flag and 
   * simplices non-split neighbors are not updated.
   * \param s    the simplex to split
   * \param seg  the index of the edge of \a s to split
   * \param[out] newVertex An array of NSEG new vertices, one for each split segment
   * \param[out] splitSimplices The two simplices generated by spliting \a s. 
   * More specifically, \a splitSimplices[0] contains the 
   * simplex generated by shrinking \a s while \a splitSimplices[1] contains 
   * its newly created partner. 
   * \param updateVertexData Whether to update the vertex data (may save time if updated
   * data is not needed).
   * \param updateSimplexData Whether to update the simplex data (may save time if updated
   * data is not needed).   
   * \tparam UpdateDataIndex If UpdateDataIndex>=0 and updateSimplexData=true, then 
   * update only the simplex data field with index 'UpdateDataIndex'. Can be used
   * to save to time ...
   * \warning Vertices' boundary flag and simplices non-split 
   * neighbors are not updated!
   */
  template <int UpdateDataIndex=-1>
    void simulateIsolatedSimplexSplit(Simplex *s, int seg,
				      Vertex &newVertex,
				      Simplex splitSimplices[2],
				      bool updateVertexData=true,
				      bool updateSimplexData=true) const
  {   
    const long rBufferSize=Vertex::template 
      getSimplexRefineBufferSize<MyType,SegmentHandle>();

    SegmentHandle handle=s->getSegmentHandle(seg);  
    char buffer[rBufferSize];

    // Create the new vertex
    if (updateVertexData)
      {
	// Update coordinates and data
	Vertex::Data::onRefineVertex(this,handle,&newVertex,&s,1,buffer);
      }
    else
      {
	// Only the coordinates of the vertex are updated in that case
	newVertex.setRefinedCoords(this,handle,&s,1,buffer);
      }	  
    
    splitSimplices[0]=(*s);    

    char *curBuffer = buffer;
    
    splitSimplices[0].splitSegment
      (handle,NULL,&newVertex,&splitSimplices[1]);

    // update simplex data if required
    if (updateSimplexData)
      {
	Simplex *s0=&splitSimplices[0];
	Simplex *s1=&splitSimplices[1];	   	      
	      
	if (UpdateDataIndex>=0)
	  {
	    // Update only data with index UpdateDataIndex
	    Simplex::Data::template refineSimplexData<UpdateDataIndex>
	      (this,handle,&newVertex,&s0,&s1,1,curBuffer);
	  }
	else
	  {
	    // update all data
	    Simplex::Data::onRefineSimplices
	      (this,handle,&newVertex,&s0,&s1,1,curBuffer);
	  }
      }
  }

  /** \brief Unrefine previously split simplices - cf. refine() -
   *
   * This function may only merge simplices that are the result of a split from refine !
   * Note that as a positive side effect, after calling coarsen the ghost/shadow 
   * layer is reduced to a minimum
   * \param solver a pointer to an object implementing a function 
   * \verbatim
     bool solver->checkCoarsen(Simplex *s1, Simplex *s2, Vertex *v1, Vertex *v2, Vertex *vSplit)
     \endverbatim
   * that returns true if simplices \a s1 and \a s2 should be merged by deleting 
   * vertex \a vSplit. For convenience, vertices \a v1 and \a v2 are the vertices that will
   * be joined by a new edge after vSplit is removed (i.e. vSplit was the vertex used to 
   * split segment [v1,v2] in function refine ). Note that this function will be called for
   * every pairs of currently existing simplices {s1,s2} that are the result of a 
   * previous split.
   * \param solver a pointer to an object implementing a function 
   * checkCoarsen. See description hereabove for more information.
   */
  template <class C>
  long coarsen(C *solver)
  {      
    if (glb::debug>1)
      dumpToNDnetwork("before_coarsening",
		      IO::NDNET_WithShadows|
		      IO::NDNET_WithGhosts|
		      IO::NDNET_WithNeighbors,
		      true);
      
    //int pass=1;
    long nCoarsened=0;
    long nCoarsenedShared=0;
    bool ghostLayerUpdated=false;
    bool shadowLayerUpdated=false;
    const int myRank=mpiCom->rank();
    const int nParts = mpiCom->size();

    if ((!glb::console->willPrint<LOG_INFO>())&&(glb::console->willPrint<LOG_STD>())) 
      glb::console->printFlush<LOG_STD>("Coarsening ... "); 

    //glb::console->printFlush<LOG_INFO>("Coarsening the mesh (P%d) ... ",pass++);
    glb::console->printFlush<LOG_INFO>("Coarsening ... ");

    glb::console->printFlush<LOG_PEDANTIC>("(candidates) ");

    // FIXME: use static to spare reallocations ?
    DICE_MESH_STATIC std::vector< Simplex* > candidates;
    const int nRequestsTotal=9;
    std::vector<MPI_Request> allRequests[nRequestsTotal];
  
    // Candidates are static, so clean-up results from previous call
    // Note that this should not deallocate the memory, which is good :)
    candidates.resize(0); 
    if (candidates.capacity()==0)
      candidates.reserve(LocalMesh::getNVertices()/10);
                
    // We need the vertices flags cleaned up !
    LocalMesh::resetVerticesTags();
   
    // First, for each simplex S, we tag any of its vertices that cannot possibly be 
    // deleted without being in conflict with the fact that S can only be merged with 
    // its partner by removing a given vertex.
    // FIXME invert the tagging to lower the computational time ?
    // FIXME writing the tags does not make openMP happy ;)
    //FIXME #pragma omp parallel for
    //for (long th=0;th<glb::num_omp_threads;++th)
      {	
	const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=it_end;++it)
     // for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,glb::num_omp_threads);
     //      it!=it_end;++it)
	  {
	    Simplex *cur=*it;
	    Simplex *partner=cur->getPartnerNotShared();
	    
	    if (partner == NULL)
	      {
		// This simplex is unrefined or its partner has been further refined, 
		// so we cannot delete ANY of its vertices
		for (int i=0;i<Simplex::NVERT;++i) 
		  cur->getVertex(i)->setTaggedF();
		// cur->cache.ptr = NULL;
	      }
	    else if (partner>cur) // so that we don't inspect the simplices twice
	      {		
		Vertex *vRef=cur->getVertex(cur->getSplitIndex());

		// cur->cache.ptr and partner->cache.ptr store a pointer to the segment
		// that would be created if cur and partner were to be merged
		for (int i=0;i<Simplex::NVERT;++i)
		  {
		    Vertex *v=cur->getVertex(i);
		    if (v!=vRef) 
		      {
			v->setTaggedF();
			if (cur->getNeighbor(i)==partner) 
			  cur->cache.ptr = v;
		      }
		  }
		for (int i=0;i<Simplex::NVERT;++i)
		  {
		    Vertex *v=partner->getVertex(i);
		    if (v!=vRef) 
		      {
			v->setTaggedF();
			if (partner->getNeighbor(i)==cur) 
			  partner->cache.ptr = v;
		      }
		  }
	      }		     
	  }
      }
  
    // Check wether the local shared vertices are flagged consistently on every process
    if (nParts>1)
      {
	glb::console->printFlush<LOG_PEDANTIC>("(com) ");
	
	//MPI_Request req[ghostExchange.sendRank.size()];
	allRequests[0].resize(ghostExchange.sendRank.size());
	std::vector< std::vector<unsigned char> > 
	  sendCancelledSplit(ghostExchange.sendRank.size());

	#pragma omp parallel for
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  {
	    const std::vector<Simplex *> &curSend = 
	      ghostExchange.send[ghostExchange.sendRank[i]];

	    sendCancelledSplit[i].resize(curSend.size());
	    // send any removable but cancelled vertex that is local and shared.
	    for (unsigned long j=0;j<curSend.size();++j)
	      {	
		unsigned char result=0;	       
		for (int k=0;k<Simplex::NVERT;++k)
		  {
		    Vertex *v=curSend[j]->getVertex(k);
		    // Vertex is shared as it belongs to a ghost simplex, just ensure that
		    // it is local
		    //if (!v->isShadowOrGhost()) 
		    if (v->isShared())
		      {
			int spi=curSend[j]->getSplitIndex();
			// if it is tagged, send the info ...
			if (k==spi)
			  {
			    if (v->isTagged()) result |= (1<<k);
			  }
			else result |= (1<<k);
		      }
		  }
		sendCancelledSplit[i][j]=result;
	      }
	    mpiCom->Isend(&sendCancelledSplit[i][0],sendCancelledSplit[i].size(),
	    		  ghostExchange.sendRank[i],&allRequests[0][i],
			  mpiTagsStart_coarsen+0);
	  }

	if (glb::debug>1)
	  dumpToNDnetwork("during_coarsening",
			  IO::NDNET_WithShadows|
			  IO::NDNET_WithGhosts|
			  IO::NDNET_WithNeighbors,
			  true);

	
	// receive everything and cancel appropriate candidates	
	// we also reset the ghosts simplices cache here to be sure we can use it later
	DICE_MESH_STATIC std::vector<unsigned char> received;
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {	
	    int count=0;
	    int source=mpiCom->ProbeCount<unsigned char>(count,mpiTagsStart_coarsen+0);
	    received.resize(0); // avoid unnecessary copies
	    received.resize(count);	
	    mpiCom->Recv(&received[0],received.size(),source,mpiTagsStart_coarsen+0);
	    const unsigned long N = received.size();
	    for (unsigned long j=0;j<N;++j)
	      {	
		for (int k=0;k<Simplex::NVERT;++k)
		  {
		    if (received[j]&(1<<k))
		      ghostExchange.receive[source][j]->getVertex(k)->setTaggedF();
		  }
		
		// We are already here, so let's reset the ghost simplices cache for later
		ghostExchange.receive[source][j]->cache.ptr=NULL;
	      }
	  }
	glb::console->printFlush<LOG_PEDANTIC>("(candidates) ");
      }

    
    // Now we want to create a list of simplices for each vertex that can be coarsened
    // without conflict. 
    // N.B.: For each vertex, we only want ONE simplex.
#pragma omp parallel for
    for (long th=0;th<glb::num_omp_threads;++th)
      {	
	const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,glb::num_omp_threads);
	     it!=it_end;++it)
	  {
	    Simplex *cur=*it;
	    Simplex *partner=cur->getPartnerNotShared();

	    if ((partner!=NULL)&&(partner>cur))
	      {	
		// check if the vertex it would refine is available
		Vertex *vRef=cur->getVertex(cur->getSplitIndex());
#pragma omp critical
		if (!vRef->isTagged()) 
		  {
		    // not tagged => we have a candidate for coarsening
		    candidates.push_back(cur);
		    // tag vRef so that we have only one simplex associated to it
		    vRef->setTaggedF();
		  }	
	      }
	  }
      }

  
    // These will store the local ids of the deleted vertices and simplices
    // This will be usefull for updating the indices
    DICE_MESH_STATIC std::vector<LocalIndex> removedVerticesLID;
    DICE_MESH_STATIC std::vector<LocalIndex> removedSimplicesLID;
    
    // This is used to store pointers to shared vertices that got deleted,
    // as we won't be able to query them anymore ...
    DICE_MESH_STATIC typename my_dense_set<Vertex*>::type rmSharedVerticesSet;
    DICE_MESH_STATIC bool rmSharedVerticesSet_init=true;
    if (rmSharedVerticesSet_init)
      {
	set_set_empty_key(rmSharedVerticesSet,NULL);
	rmSharedVerticesSet_init=false;
      }
    rmSharedVerticesSet.clear_no_resize();
   

    // candidates are static, so clean-up results from previous call
    // Note that this should not deallocate the memory, which is good :)
    removedSimplicesLID.resize(0);
    removedVerticesLID.resize(0);

    // Try to estimate the sizes of the array if it is too small
    if (removedVerticesLID.capacity() < 2*candidates.size())     
      {
	removedVerticesLID.reserve(candidates.size());
	if (NDIM==1) removedSimplicesLID.reserve(candidates.size());
	if (NDIM==2) removedSimplicesLID.reserve(candidates.size()*2);
	else removedSimplicesLID.reserve(candidates.size()*5); // wild guess ;)
      }

    //candidates.resize(0);

    typedef typename my_dense_hash<Simplex*,Simplex*>::type SimplexSimplexHash;
    typedef typename SimplexSimplexHash::iterator SimplexSimplexHash_it;
    DICE_MESH_STATIC SimplexSimplexHash modifiedSimplicesHash;
    DICE_MESH_STATIC bool simplexHashInit=true;
    if (simplexHashInit)
      {
	set_hash_empty_key(modifiedSimplicesHash,static_cast<Simplex*>(NULL));
	modifiedSimplicesHash.rehash(LocalMesh::getNGhostSimplices());
	simplexHashInit=false;
      }
    modifiedSimplicesHash.clear_no_resize();
    
    if (candidates.size()>0)
      {
	glb::console->printFlush<LOG_PEDANTIC>("(removing) ");		

	// Now that we have the list of all the vertices that could be deleted, check
	// whether we should remove them.
	// FIXME: use openMP here ?
#pragma omp parallel for
	for (long th=0;th<glb::num_omp_threads;++th)
	  {
	    const unsigned long delta=candidates.size()/glb::num_omp_threads;
	    const unsigned long i0=th*delta;
	    const unsigned long imax=
	      (th==(glb::num_omp_threads-1))?candidates.size():i0+delta;
	   
	    // This will be used to store all the simplices that will be merged when 
	    // removing a vertex. We declare them here to avoid reallocations.
	    std::vector<Simplex*> s1;
	    std::vector<Simplex*> s2;
	    for (unsigned long i=i0;i<imax;++i)
	      {			     
		// Two partners sharing the vertex to be removed
		Simplex *cur=candidates[i];		
		Simplex *partner=cur->getPartnerNotShared();
		// the vertex to be removed
		Vertex *v=cur->getVertex(cur->getSplitIndex());

		// [curV,partnerV] would be the segment created if v was deleted
		Vertex *curV=static_cast<Vertex*>(cur->cache.ptr);	   
		Vertex *partnerV=static_cast<Vertex*>(partner->cache.ptr);
		
		// 'dirV' is used to ensure we are turning in the same direction around
		// both segments.
		Vertex *dirV;
		for (int j=0;j<Simplex::NVERT;++j)
		  {
		    dirV = cur->getVertex(j);
		    if ((dirV!=v)&&(dirV!=curV))
		      break;
		  }
		
		// Circulators around the merging segments
		segment_circulator ci_end=Segment(v,curV,cur).getCirculator(dirV);
		segment_circulator ci=ci_end;
		segment_circulator pci_end=Segment(v,partnerV,partner).getCirculator(dirV);
		segment_circulator pci=pci_end;
	    
		s1.resize(0);
		s2.resize(0);
		// lists of simplices affected by the removal of v
		// If v is deleted, then s1[i] and s2[i] should merge		
		do
		  {
		    s1.push_back(*ci);
		    s2.push_back(*pci);			   
		    ++pci;
		  }
		while ((++ci)!=ci_end);
		
		// Now we can check wether the solver wants to coarsen ...
		if (solver->checkCoarsen(s1,s2,curV,partnerV,v))
		  {
		    // YES -> remove this vertex by merging s1[n] with s2[n]
#pragma omp critical
		    {		    		    		   		    
		      for (unsigned long n=0;n<s2.size();++n)
			{		
			  Simplex *rm=s2[n];
			  Simplex *keep=s1[n];	

			  if (s2[n]->isShadowOrGhost())
			    continue;		   			  

			  // Select which one is deleted, which one is expanded.
			  // We let Simplex::merge decide so that it is done consistently
			  rm = s2[n]->merge(Tree::nodePool);
			  keep = (s1[n]==rm)?s2[n]:s1[n];
			  int keepIndex=-1;
			  // update the neighbors' neighbors
			  for (int j=0;j<Simplex::NNEI;++j)
			    {
			      Simplex *nei=rm->getNeighbor(j);
			      if (nei!=NULL)
				{
				  if (nei!=keep)
				    {
				      int k=nei->getNeighborIndex(rm);
				      // neighbors's neighbors may already have been 
				      // updated so we have to check that 'rm' is still 
				      // considered a neighbor by its neighbor ...
				      if (k>=0) nei->setNeighbor(k,keep);
				    }
				  else keepIndex = j;
				}
			    }

			  // update 'keep' neighbors
			  int rmIndex=keep->getNeighborIndex(rm);
			  int rmVId = rm->getVertexIndex(v);
			  Simplex* newNei=rm->getNeighbor(rmVId);
			  keep->setNeighbor(rmIndex,newNei);
		  
			  // Update keep's simplices data
			  int keepVId = keep->getVertexIndex(v);
			  
			  Simplex::Data::template onCoarsenSimplex<MyType,Simplex,Vertex>
			    (this,keep,keepVId,rmIndex,rm,rmVId,keepIndex,
			     static_cast<Vertex*>(rm->cache.ptr));

			  // solver->onCoarsen(keep,keepVId,rmIndex,
			  // 		    rm,rmVId,keepIndex,
			  // 		    static_cast<Vertex*>(rm->cache.ptr));

			  // and replace the deleted vertex in 'keep'
			  keep->setVertex(keepVId,static_cast<Vertex*>(rm->cache.ptr));
			  keep->setGeneration(keep->getGeneration().rank()-1,
					      keep->getGeneration().id());

			  if (nParts>1)
			    {
			      // we'll need that info to communicate with other processes
			      keep->cache.c[0]=keepVId; // removed vertex
			      keep->cache.c[1]=rmIndex; // partner
			      keep->cache.c[2]=rmVId; // removed vertex in rm
			      keep->cache.c[3]=keepIndex; //partner in rm 
			      modifiedSimplicesHash.insert(std::make_pair(keep,keep));
			      modifiedSimplicesHash.insert(std::make_pair(rm,keep));
			    }
			  
			  removedSimplicesLID.push_back(rm->getLocalIndex());
			  LocalMesh::simplexPool.recycle(rm);
			  
			}
		      // all simplices are now merged, we can delete the vertex.
		      removedVerticesLID.push_back(v->getLocalIndex());	
		      // if (glb::debug)
		      // 	v->template print<LOG_STD_ALL>();
		      if (v->isShared()) rmSharedVerticesSet.insert(v);
		      LocalMesh::vertexPool.recycle(v);			      
		    } // end critical		   
		  } // end coarsening
	      } // end check candidates	    
	  } // end openMP loop

      } // candidates.size() > 0
     
    // Declare the structures used to ISend : because we are transfering in the 
    // background, we would not want them to be deallocated before finishing, so to avoid
    // having to synchronize in the 'if (nparts>1)' section, we declare them here and 
    // will only wait for them to complete at the end of coarsen()
   
    std::vector< std::vector<MpiExchg_CoarsenGhosts> > 
      sentGhosts(ghostExchange.receiveRank.size());       
    // std::vector< std::vector<MpiExchg_CoarsenShadows> > 
    //   sentShadows(shadowExchange.receiveRank.size()); 
    std::vector< std::vector<char> > 
      sentGhostsUpdate(ghostExchange.receiveRank.size());
    std::vector< std::vector<MpiExchg_NewShadowSimplices> >
      sentNewShadowSimplices(ghostExchange.sendRank.size());
    std::vector< std::vector<GlobalIdentityValue> > 
      queryTransfer(ghostExchange.sendRank.size());

    // we may need these hash tables to retrieve the shadow vertices/simplices 
    // that we will have to create
    typedef typename my_dense_hash<GlobalIdentityValue,Simplex*>::type 
      GidShadowSimplexHash;
    typedef typename my_dense_hash<GlobalIdentityValue,Vertex*>::type 
      GidShadowVertexHash;
    DICE_MESH_STATIC GidShadowSimplexHash shadowSimplexHash;
    DICE_MESH_STATIC GidShadowVertexHash shadowVertexHash;

    if (nParts>1)
      {
	glb::console->printFlush<LOG_PEDANTIC>("(com) ");  

	std::vector< std::vector<MpiExchg_CoarsenGhosts> > 
	  receivedGhosts(ghostExchange.receiveRank.size());

	//FIXME : optimize with irecv/iSend ?
	// for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	//   {
	//     receivedGhosts[source].resize(ghostExchange.receive[i].size());
	//     mpiCom->IRecvMpiType(&receivedGhosts[i].front(),receivedGhosts[i].size(),
	// 			 mpiType_coarsenGhosts,ghostExchange.receiveRank[i],1);
	//   }
       
	allRequests[1].resize(ghostExchange.sendRank.size());
#pragma omp parallel for
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  {
	    const std::vector<Simplex *> &curSend = 
	      ghostExchange.send[ghostExchange.sendRank[i]];
	    MpiExchg_CoarsenGhosts val;
	    sentGhosts[i].reserve(curSend.size()/2);
	    for (unsigned long j=0;j<curSend.size();++j)
	      {
		SimplexSimplexHash_it it = modifiedSimplicesHash.find(curSend[j]);
		if (it != modifiedSimplicesHash.end())
		  {
		    //if (LocalMesh::simplexPool.isFree(curSend[j]))
		    if (it->first!=it->second)
		      {
			// This is the deleted simplex in the pair
			// We send it in case its partner is a shadow, as we would otherwise
			// miss it ...
			val.set(it->second,j,0);		
			sentGhosts[i].push_back(val);			 
		      }
		    else
		      {
			// this is the coarsened simplex in the pair
			val.set(curSend[j],j,1);
			sentGhosts[i].push_back(val);
		      }			    
		  }
	      }
	    // glb::console->print<LOG_STD_ALL>("process %d sending to %d.\n",
	    // 				     mpiCom->rank(),ghostExchange.sendRank[i]);
	    mpiCom->IsendMpiType(&sentGhosts[i][0],
				 sentGhosts[i].size(),
				 mpiType_coarsenGhosts.getType(),
				 ghostExchange.sendRank[i],
				 &allRequests[1][i],mpiTagsStart_coarsen+1);
	  } 

	// for (unsigned long i=0;i<shadowExchange.receiveRank.size();++i)
	//   {
	//     receivedShadows[source].resize(shadowExchange.receive[i].size());
	//     mpiCom->IRecvMpiType(&receivedShadows[i].front(),receivedShadows[i].size(),
	// 			 mpiType_coarsenShadows,shadowExchange.receiveRank[i],2);
	//   }
	/*
#pragma omp parallel for
	for (unsigned long i=0;i<shadowExchange.sendRank.size();++i)
	  {
	    const std::vector<Simplex *> &curSend = 
	      shadowExchange.send[shadowExchange.sendRank[i]];
	    MpiExchg_CoarsenShadows val;
	    sentShadows[i].reserve(curSend.size()/2);
	    for (unsigned long j=0;j<curSend.size();++j)
	      {
		SimplexSimplexHash_it it = modifiedSimplicesHash.find(curSend[j]);
		if (it != modifiedSimplicesHash.end())
		  {
		    //if (LocalMesh::simplexPool.isFree(curSend[j]))
		    if (it->first!=it->second)
		      {
			// This is the deleted simplex in the pair
			val.set(it->second,j,0);		
			sentShadows[i].push_back(val);
		      }
		    else
		      {
			// this is the coarsened simplex in the pair
			val.set(curSend[j],j,1);
			sentShadows[i].push_back(val);		
		      }	
		  }
	      }	   
	    mpiCom->IsendMpiType(&sentShadows[i][0],
				 sentShadows[i].size(),
				 mpiType_coarsenShadows,
				 shadowExchange.sendRank[i],
				 &sentShadowsRequests[i],2);
	  } 
	*/


	std::map<int,int> receiveGhostRankIndex;    
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  receiveGhostRankIndex.insert
	    (std::make_pair(ghostExchange.receiveRank[i],i));
	//receiveGhostRankIndex[ghostExchange.receiveRank[i]] = i;
	std::map<int,int> receiveShadowRankIndex;    
	for (unsigned long i=0;i<shadowExchange.receiveRank.size();++i)
	  receiveShadowRankIndex.insert
	    (std::make_pair(shadowExchange.receiveRank[i],i));

	std::map<int,int> sendGhostRankIndex;    
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  sendGhostRankIndex.insert
	    (std::make_pair(ghostExchange.sendRank[i],i));
	std::map<int,int> sendShadowRankIndex;    
	for (unsigned long i=0;i<shadowExchange.sendRank.size();++i)
	  sendShadowRankIndex.insert
	    (std::make_pair(shadowExchange.sendRank[i],i));

	// FIXME: a simple map should be good enough I guess ?
	std::map<ShadowVertex*,GhostVertex*> ShadowVertex2GhostVertexMap;

	// FIXME: don t need the receivedGhosts here ? => use local vector
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    int count=0;
	    int source=mpiCom->ProbeCount(count,mpiType_coarsenGhosts.getType(),
					  mpiTagsStart_coarsen+1);
	    int srcIndex=receiveGhostRankIndex.at(source);
	    
	    if (count>0) receivedGhosts[srcIndex].resize(count);
	    else receivedGhosts[srcIndex].clear();	    

	    mpiCom->RecvMpiType(&receivedGhosts[srcIndex][0],count,
				mpiType_coarsenGhosts.getType(),source,
				mpiTagsStart_coarsen+1);
	    
	    for (unsigned long j=0;j<receivedGhosts[srcIndex].size();++j)
	      {
		unsigned int index = receivedGhosts[srcIndex][j].getBaseCellIndex();
		
		if (receivedGhosts[srcIndex][j].type==0) 
		  {	
		    // This simplex is the deleted one in the pair, but its partner 
		    // may be a shadow in which case the result would be a shadow instead
		    // of a ghost.
		    GhostSimplex *rm_= ghostExchange.receive[source][index];
		    Simplex *keep_=rm_->getNeighbor(receivedGhosts[srcIndex][j].
						    partnerIndex);
		    if (!keep_->isShadow()) continue;

		    // NB: we will first swap keep and rm so that rm is actually coarsened
		    // and transformed into a GhostSimplex while keep will be deleted
		    // Note that we MUST keep the vertices and neighbors order as though
		    // we had not swapped.
		    Vertex *v=rm_->getVertex(receivedGhosts[srcIndex][j].rmVertexIndex);
		    GhostSimplex *keep = rm_;
		    ShadowSimplex *rm = static_cast<ShadowSimplex*>(keep_);

		    GlobalIdentityValue tmp= keep->getGlobalIdentity().get();
		    keep->setGlobalIdentity(rm->getGlobalIdentity());
		    keep->setData(receivedGhosts[srcIndex][j].sData);
		    rm->setGlobalIdentity(tmp);
		    
		    // check the order of vertices that 'keep' should have.
		    // By changing the order that way, we have the same result as
		    // if rm had been coarsened and keep discarded !
		    int vertOrder[Simplex::NVERT];		    		    
		    for (int i=0;i<Simplex::NVERT;++i)
		      {
			if (rm->getVertex(i) == v)
			  vertOrder[i]=keep->getNeighborIndex(rm);
			else
			  vertOrder[i]=keep->getVertexIndex(rm->getVertex(i));
			if (vertOrder[i]<0)
			  vertOrder[i]=keep->getVertexIndex(v);
		      }
		    
		    //v->template print<LOG_STD_ALL>();
		    Vertex* newV=NULL; // the vertex that will replace 'v' in keep
		    // update the neighbors' neighbors
		    for (int k=0;k<Simplex::NNEI;++k)
		      {
			Simplex *nei=rm->getNeighbor(k);
			if (nei!=NULL)
			  {
			    if (nei!=keep)
			      {			
				int l=nei->getNeighborIndex(rm);
				if (l>=0) nei->setNeighbor(l,keep);
			      }
			    else newV=rm->getVertex(k);
			  }
		      }

		     // update 'keep' neighbors
		    int newNeiId=keep->getNeighborIndex(rm);
		    Simplex* newNei=rm->getNeighbor(rm->getVertexIndex(v));
		    if (newNei != NULL)
		      keep->setNeighbor(newNeiId,newNei);
		    else
		      {
			// This may happen when rm is a shadow as we do not know
			// all the neighbors of a shadowSimplex !
			// In that case, we cannot trust the NULL status of the new 
			// neighbor and set it to 'keep' itself to tag it
			keep->setNeighbor(newNeiId,keep);
		      }
		   	
		    // and replace the deleted vertex in 'keep'	
		    Vertex *newGV=newV;
		    if (newV->isShadow())
		      {
			ShadowVertex *sv = static_cast<ShadowVertex*>(newV);
			// now v will belong to a ghost simplex, so it should be a
			// ghost vertex !
			typename std::map<ShadowVertex*,GhostVertex*>::iterator it =
			  ShadowVertex2GhostVertexMap.find(sv);
			
			if (it == ShadowVertex2GhostVertexMap.end())
			  {
			    GhostVertex *gv;
			    LocalMesh::ghostVertexPool.pop(&gv);			    
			    gv->copy(newV);
			    ShadowVertex2GhostVertexMap.insert(std::make_pair(sv,gv));
			    newGV = static_cast<Vertex*>(gv);
			  }
			else newGV = static_cast<Vertex*>(it->second);
		      }
		    keep->setVertex(keep->getVertexIndex(v),newGV);

		    // tag for removal !
		    rm->setSetF(false);		   
		    v->setSetF(false);		  
		    rm->cache.c[0]=2; // means it was swapped
		    keep->cache.c[0]=2; // means it was swapped
		   
		    Vertex *vertices[Simplex::NVERT];
		    Simplex *neighbors[Simplex::NNEI];
		    // and reorder the vertices/neighbors
		    for (int i=0;i<Simplex::NVERT;++i)
		      {
			vertices[i]=keep->getVertex(vertOrder[i]);
			neighbors[i]=keep->getNeighbor(vertOrder[i]);
		      }
		    for (int i=0;i<Simplex::NVERT;++i)
		      {
			keep->setVertex(i,vertices[i]);
			keep->setNeighbor(i,neighbors[i]);
		      }
		    // glb::console->print<LOG_STD_ALL>("OOOOOOO:");
		    // keep->template print<LOG_STD_ALL>();
		    // rm->template print<LOG_STD_ALL>();
		    nCoarsenedShared++;		    
		  }
		else 
		  {		    
		    //coarsened partner
		    GhostSimplex *keep=ghostExchange.receive[source][index];
		    int pIndex=receivedGhosts[srcIndex][j].partnerIndex;
		    int rmVId=receivedGhosts[srcIndex][j].rmVertexIndex;
		    Simplex *rm=keep->getNeighbor(pIndex);		    
		    Vertex *v=keep->getVertex(rmVId);
		    //bool follow=false;

		    keep->setData(receivedGhosts[srcIndex][j].sData);
		    
		    // if ((keep->getGlobalIdentity().rank()==3)&&
		    // 	(keep->getGlobalIdentity().id()==24))
		    //   {
		    // 	follow=true;
		    // 	keep->template print<LOG_STD_ALL>();
		    // 	rm->template print<LOG_STD_ALL>();			
		    //   }

		    Vertex* newV=NULL; // the vertex that will replace 'v' in keep
		    // update the neighbors' neighbors
		    for (int k=0;k<Simplex::NNEI;++k)
		      {
			Simplex *nei=rm->getNeighbor(k);
			if (nei!=NULL)
			  {
			    if (nei!=keep)
			      {			
				int l=nei->getNeighborIndex(rm);
				if (l>=0) nei->setNeighbor(l,keep);
			      }
			    else newV=rm->getVertex(k);
			  }
		      }

		    // update 'keep' neighbors
		    int newNeiId=keep->getNeighborIndex(rm);
		    Simplex* newNei=rm->getNeighbor(rm->getVertexIndex(v));		   
		    if ((newNei != NULL)||(!rm->isShadow()))
		      keep->setNeighbor(newNeiId,newNei);
		    else
		      {
			// This may happen when rm is a shadow as we do not know
			// all the neighbors of a shadowSimplex !
			// In that case, we cannot trust the NULL status of the new 
			// neighbor and set it to 'keep' itself to tag it
			keep->setNeighbor(newNeiId,keep);
		      }

		    Vertex *newGV=newV;
		    if (newV->isShadow())
		      {
			ShadowVertex *sv = static_cast<ShadowVertex*>(newV);
			// now v will belong to a ghost simplex, so it should be a
			// ghost vertex !
			typename std::map<ShadowVertex*,GhostVertex*>::iterator it =
			  ShadowVertex2GhostVertexMap.find(sv);
			
			if (it == ShadowVertex2GhostVertexMap.end())
			  {
			    GhostVertex *gv;
			    LocalMesh::ghostVertexPool.pop(&gv);			    
			    gv->copy(newV);
			    ShadowVertex2GhostVertexMap.insert(std::make_pair(sv,gv));
			    newGV = static_cast<Vertex*>(gv);
			  }
			else newGV = static_cast<Vertex*>(it->second);
		      }
		    //keep->setVertex(keep->getVertexIndex(v),newGV);
		    
		    // Replace the deleted vertex in 'keep'	
		    // NB: if rm->isShadow(), newV may wrongly remain a shadowVertex
		    keep->setVertex(rmVId,newGV);
		    keep->cache.c[0]=1; // means it was kept
		    rm->cache.c[0]=1; // means it was kept
		    // tag for removal !
		    rm->setSetF(false);	
		    // check if the vertex is not deleted already (this is the case for 
		    // the removed shared vertices)
		    if (rmSharedVerticesSet.find(v) == rmSharedVerticesSet.end())
		      v->setSetF(false);
		    nCoarsenedShared++;
		    // if (follow)
		    //   {
		    // 	keep->template print<LOG_STD_ALL>();
		    // 	rm->template print<LOG_STD_ALL>();
		    //   }
		  }
	      }
	  } // ghostExchange.receive


	// std::vector< std::vector<MpiExchg_CoarsenShadows> > 
	//   receivedShadows(shadowExchange.receiveRank.size());
	
	/*
	// FIXME: don t need the receivedShadows here ? => use local vector
	// FIXME: just delete anything tagged !
	for (unsigned long i=0;i<shadowExchange.receiveRank.size();++i)
	  {
	    int count=0;
	    int source=mpiCom->ProbeCount(count,mpiType_coarsenShadows,2);
	    int srcIndex=receiveShadowRankIndex[source];
	    
	    if (count>0) receivedShadows[srcIndex].resize(count);
	    else receivedShadows[srcIndex].clear();	    

	    mpiCom->RecvMpiType(&receivedShadows[srcIndex][0],count,
				mpiType_coarsenShadows,source,2);
	    
	    for (unsigned long j=0;j<receivedShadows[srcIndex].size();++j)
	      {
		unsigned int index = receivedShadows[srcIndex][j].getBaseCellIndex();
		if (receivedShadows[srcIndex][j].type==0) 
		  {
		    // this one was deleted !
		    ShadowSimplex *rm= shadowExchange.receive[source][index];
		    //Simplex *keep=rm->getNeighbor(receivedGhosts[srcIndex][j].
		    // 				  partnerIndex);
		    Vertex *v=rm->getVertex(receivedShadows[srcIndex][j].rmVertexIndex);

		    v->setSetF(false);
		    rm->setSetF(false);
		    for (int k=0;k<Simplex::NVERT;++k)
		      {
			Simplex *s=rm->getNeighbor(k);
			if (s!=NULL)
			  {
			    int id=s->getNeighborIndex(rm);
			    if (id>=0) s->setNeighbor(id,s);
			  }
		      }
		    nCoarsenedShared++;
		  }
		else
		  {
		    // this shadow is coarsened
		    ShadowSimplex *keep= shadowExchange.receive[source][index];
		    int pIndex=receivedShadows[srcIndex][j].partnerIndex;
		    int rmVId=receivedShadows[srcIndex][j].rmVertexIndex;
		    Simplex *rm=keep->getNeighbor(pIndex);		    
		    Vertex *v=keep->getVertex(rmVId);
		    
		    // if keep->isSet(), then this was a pair of ghost/shadow where 
		    // the shadow was supposed to be deleted, but we swapped instead:)
		    //if ((rm==NULL)||(keep->isSet()))
		      {
			// we do not have enough info on the neighbor to efficiently merge
			// it, so remove it so that it gets fully transfered later
			keep->setSetF(false);
			v->setSetF(false);
			for (int k=0;k<Simplex::NVERT;++k)
			  {
			    Simplex *s=keep->getNeighbor(k);
			    
			    if ((s!=NULL)&&(s!=keep))
			      {
				int id=s->getNeighborIndex(keep);
				if (id>=0) s->setNeighbor(id,s);
			      }
			  }
			nCoarsenedShared++;				
		      }	
		      	      
		    // else
		    //   {
		    // 	// luckily we already have the new neighbor locally, so just delete
		    // 	// the vertex
		    // 	v->setSetF(false);
		    //   }	
		      
		      	   
		  }
	      }
	  } //shadowExchange.receive
	*/

	// Initialize the shadowVertices hash and shadowSimplices hash if necessary, and clear
	// them without freeing memory
	DICE_MESH_STATIC bool shadowHashInit=true;
	if (shadowHashInit)
	  {
	    set_hash_empty_key(shadowSimplexHash,GlobalIdentity::empty.get());
	    set_hash_empty_key(shadowVertexHash,GlobalIdentity::empty.get());
	    set_hash_deleted_key(shadowVertexHash,GlobalIdentity::max.get());
	    shadowSimplexHash.rehash(LocalMesh::getNShadowSimplices()/
				(shadowSimplexHash.max_load_factor()-0.01));
	    shadowVertexHash.rehash(LocalMesh::getNShadowVertices()/
				(shadowVertexHash.max_load_factor()-0.01));
	    shadowHashInit=false;
	  }
	shadowSimplexHash.clear_no_resize();
	shadowVertexHash.clear_no_resize();	

	// static Simplex* noUse;
	// noUse=NULL;

	allRequests[3].resize(ghostExchange.receiveRank.size());

	// First tag for removal any unnecessary ghost simplex.
	// => An unnecessary ghost simplex does not have any local vertex anymore.
	// FIXME: openMP ?
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    std::vector<GhostSimplex *> &curSend = 
	      ghostExchange.receive[ghostExchange.receiveRank[i]];
	    sentGhostsUpdate[i].resize(curSend.size());
	  	    
	    for (unsigned long j=0;j<curSend.size();++j)
	      {
		GhostSimplex *s=curSend[j];		
		if (s->isSet())
		  {
		    // Check if this ghost really needs to remain a ghost (i.e. it still
		    // has at least one local vertex)
		    bool remove=true;
		    for (int k=0;k<Simplex::NVERT;++k)
		      if (!s->getVertex(k)->isShadowOrGhost()) 
			remove=false;
		  
		    if (remove) 
		      s->setSetF(false);
		    else
		      {
			// Prepare the shadow vertex map, as some shadow vertex may already
			// be local as ghosts
			for (int k=0;k<Simplex::NVERT;++k)
			  {
			    Vertex *v=s->getVertex(k);
			    if (v->isGhost())
			      shadowVertexHash[v->getGlobalIdentity().get()]=v;
			  } 
		      }
		  }
	      }
	  }

	// Now we want to send back information about which ghost simplices need an update
	// For instance, it can be the need for a new neighbor shadow, an update of the
	// ghostExchange structure, or a removal 
	// Note: we also remove any deleted simplex from the ghostExchange.receive structure
	// FIXME: USE OPENMP ?	
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    std::vector<GhostSimplex *> &curSend = 
	      ghostExchange.receive[ghostExchange.receiveRank[i]];
	    sentGhostsUpdate[i].resize(curSend.size());

	    unsigned long newIndex=0;
	    for (unsigned long j=0;j<curSend.size();++j)
	      {
		GhostSimplex *s=curSend[j];
		char status=0;
		if (s->isSet())
		  {
		    // This happens when a ghost was merged with a shadow, 
		    // and the ghost was supposed to be deleted.
		    if (s->cache.c[0]==2)
		      status |= (1<<6);
		    		    
		    // Check if the simplex has all the neighbors it needs
		    // We have set invalid neighbors to the simplex itself earlier
		    //Vertex *v;
		    int count=0;
		    for (int k=0;k<Simplex::NNEI;++k)
		      {
			Simplex *nei = s->getNeighbor(k);
			
			if (nei!=NULL)
			  {
			    // we will need to retrieve this neighbor full info
			    // as it will be a ShadowSimplex ...
			    //nei->cache.c[1]=1; 			    
			    if ((nei==s)||(!nei->isSet())||(nei->isShadow()))
			      {
				status |= (1<<k);
				curSend[j]->cache.c[1]=k;
				count++;
			      }
			    /*
			    else if (!nei->isSet()) 
			      {
				if (nei->isGhost()) 
				  {
				    printf("SHOULD NOT HAPPEN YET\n");exit(0);
				    // that's a deleted ghost that still exists
				    // we just need it as a shadow !
				    // GidShadowSimplexHash_it it = 
				    //   newShadowHash.find(nei->getGlobalIndex(myRank).get());
				    // if (it == newShadowHash.end())
				    //   {
				    // 	ShadowVertex *shadow;
				    // 	LocalMesh::ShadowPool.pop(&shadow);
				    // 	shadow->copy(static_cast<Ghost*>(nei));
				    // 	newShadowHash.insert(std::make_pair(nei->getGlobalIndex(myRank).get(),shadow));
				    // 	s->setNeighbor(shadow);
				    //   }
				    // else s->setNeighbor(it->second);
				  }
				else 
				  {
				    printf("SHOULD NOT HAPPEN YET\n");exit(0);
				    // this is the result of a remote merger !
				    status |= (1<<k);
				    count++;
				  }
			      }
			    */			    
			  }
		      }	
		    // if (noUse != NULL)
		    //   {
		    // 	curSend[j]->template print<LOG_WARNING>();
		    // 	noUse->template print<LOG_WARNING>();
		    //   }
		    // if (curSend[j]->getGlobalIdentity().id() == 2816)
		    //   {
		    // 	noUse=curSend[j];
		    // 	glb::console->print<LOG_WARNING>("Found 2816 with status %d\n",status);
		    // 	curSend[j]->template print<LOG_WARNING>();
		    //   }
		    curSend[newIndex++]=curSend[j];
		    // So we can retrieve the index of any ghost in ghostExchange.receive
		    curSend[newIndex-1]->cache.ui[1] = newIndex-1;
		    if (count>1) 
		      {
			PRINT_SRC_INFO(LOG_ERROR);
			glb::console->print<LOG_ERROR>("Found more than one unresolved neighbor for a ghost simplex.\n");
			glb::console->print<LOG_ERROR>("This is never supposed to happen ...");
			exit(-1);			
		      }
		  }
		else status |= (1<<7); // Deleted

		sentGhostsUpdate[i][j]=status;
		//if (s->isSet()) curSend[newIndex++]=curSend[j];
	      }
	    if (newIndex != curSend.size()) 
	      {
		ghostLayerUpdated=true;				
		curSend.resize(newIndex);
	      }
	    
	    mpiCom->Isend(&sentGhostsUpdate[i][0],sentGhostsUpdate[i].size(),
			  ghostExchange.receiveRank[i],
			  &allRequests[3][i],mpiTagsStart_coarsen+3);
	  }

	// Now we can free all the shadows as we will rebuild them entirely
	for (long i=0;i<shadowExchange.receiveRank.size();++i)
	  shadowExchange.receive[shadowExchange.receiveRank[i]].clear();
	for (long i=0;i<shadowExchange.sendRank.size();++i)
	  shadowExchange.send[shadowExchange.sendRank[i]].clear();
	const shadowSimplexPtr_iterator itss_end=LocalMesh::shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();
	     it!=itss_end;++it)
	  LocalMesh::shadowSimplexPool.recycle(*it);
	const shadowVertexPtr_iterator itsv_end=LocalMesh::shadowVertexEnd();
	for (shadowVertexPtr_iterator it=LocalMesh::shadowVertexBegin();
	     it!=itsv_end;++it)
	  LocalMesh::shadowVertexPool.recycle(*it);
	shadowLayerUpdated=true;

	
	// Check which shadow simplices are not needed anymore
	std::vector< std::vector<char> > 
	  receivedGhostsUpdate(ghostExchange.sendRank.size());
	// std::vector< std::vector<MpiExchg_NewShadowSimplices> >
	//   sentNewShadowSimplices(ghostExchange.sendRank.size());
	std::vector< std::vector<MpiExchg_NewShadowSimplices> >
	  receivedNewShadowSimplices(ghostExchange.receiveRank.size());

	// std::vector<int> shadowSimplicesNeedUpdate(shadowExchange.receive.size(),0);
	// const shadowSimplexPtr_iterator itss_end=LocalMesh::shadowSimplexEnd();
	// for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();
	//      it!=itss_end;++it)
	//   {
	//     if (it->cache.c[1]!=1)
	//       shadowSimplicesNeedUpdate[it->getGlobalIdentity().rank()]=1;
	//   }

	// // And recycle them ... If one shadow from a given process has become
	// // useless, we recycle all the shadows from that process !
	// for (long i=0;i<shadowExchange.receiveRank.size();++i)
	//   {
	//     if (!shadowSimplicesNeedUpdate[shadowExchange.receiveRank[i]])
	//       continue;
	//     std::vector<ShadowSimplex *> &curReceive = 
	//       shadowExchange.receive[shadowExchange.receiveRank[i]];
	    
	//     for (long j=0;j<curReceive.size();++j)
	//       if (!LocalMesh::shadowSimplexPool.isFree(curReceive[j]))
	// 	LocalMesh::shadowSimplexPool.recycle(&curReceive[j]);
	//     curReceive.resize(0);
	//   }

	// std::vector< std::vector<GlobalIdentityValue> > 
	//   queryTransfer(ghostExchange.sendRank.size());

	allRequests[4].resize(ghostExchange.sendRank.size());
	// Update the ghostexchange.send structure
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  {
	    int source=mpiCom->Probe(mpiTagsStart_coarsen+3);
	    int srcIndex = sendGhostRankIndex.at(source);

	    std::vector<Simplex *> &curReceive = 
	      ghostExchange.send[ghostExchange.sendRank[srcIndex]];
	    receivedGhostsUpdate[srcIndex].resize(curReceive.size());

	    mpiCom->Recv(&receivedGhostsUpdate[srcIndex][0],
			 curReceive.size(),source,mpiTagsStart_coarsen+3);
	    
	    unsigned long newIndex=0;
	    for (unsigned long j=0;j<curReceive.size();++j)
	      {
		const char status=receivedGhostsUpdate[srcIndex][j];
		if (status!=0)
		  {
		    if (status&(1<<7)) // remove it
		      {
			
		      }
		    else 
		      {
			if (status&(1<<6)) 
			  {
			    // the deleted ghost and the coarsened shadow simplices were
			    // swapped remotely, so swap it in ghostExchange.send also
			    SimplexSimplexHash_it it = 
			      modifiedSimplicesHash.find(curReceive[j]);
			    curReceive[j]=it->second;
			  }

			Simplex *nei;
			MpiExchg_NewShadowSimplices val;
			// check if we need to send back neighbor info
		     	for (int k=0;k<Simplex::NNEI;++k)
			  if (status&(1<<k))
			    {
			      // YES -> send the full simplex info, it will be a shadow.
			      // NB1: we have to be carefull here as the shadow may 
			      // originate from a different process (->transfer query)!
			      // NB2: we use newIndex, NOT j, as the ghostExchange
			      // structure is being resized !
			      val.set(curReceive[j],k,myRank,newIndex);
			      nei=curReceive[j]->getNeighbor(k);
			    }

			if (!val.isEmpty())
			  {			   
			    sentNewShadowSimplices[srcIndex].push_back(val);
			    // if ((glb::debug)&&(nei!=NULL))
			    //   {
			    // 	glb::console->print<LOG_STD_ALL>("Sending query for baseIndex %d\n",newIndex);
			    // 	curReceive[j]->template print<LOG_STD_ALL>("CUR");
			    // 	nei->template print<LOG_STD_ALL>("NEI");
			    //   }
			    // update shadowExchange.send if necessary
			    // FIXME: write this IF correctly !
			    if ((nei!=NULL)&&(!nei->isShadowOrGhost()))
			      shadowExchange.
				send[ghostExchange.sendRank[srcIndex]].
				push_back(nei);
			    else if (nei!=NULL)
			      {	
				// This simplex actually belongs to a neighbor process,
				// we have to transfer the query so that the 
				// shadowExchange structure of that process can be updated
				int dest=nei->getGlobalIdentity(myRank).rank();
				int destIndex = receiveGhostRankIndex.at(dest);
				GlobalIdentity gid(val.simplexGid);

				// Set the gid rank to the rank of the process where the
				// ghostSimplex that needs info on its shadow neighbor is
				// i.e.: the process from which we received this query
				gid.setRank(source);//ghostExchange.sendRank[srcIndex]);
				queryTransfer[destIndex].push_back(gid.get());
			
				// -> process with rank 'dest' stores this simplex in 
				// ghostExchange.receive[myRank][nei->cache.ui[1]], so it is
				// easy to retrieve it if we also send that piece of info
				queryTransfer[destIndex].push_back(nei->cache.ui[1]);
			
				// glb::console->
				//   print<LOG_STD_ALL>("QUERY TRANSFERED (%d->%d->%d):\n",
				// 		     source,myRank,dest);
				
				// curReceive[j]->template print<LOG_STD_ALL>();
				// nei->template print<LOG_STD_ALL>();
				// for (int kk=0;kk<Simplex::NVERT;++kk)
				//   {
				//     curReceive[j]->getVertex(kk)->
				//       template print<LOG_STD_ALL>();
				//   }
			      }
			  }
			curReceive[newIndex++]=curReceive[j];
		      }
		  }
		else curReceive[newIndex++]=curReceive[j];
	      }
	    if (curReceive.size() != newIndex)
	      curReceive.resize(newIndex);
	      	    
	    mpiCom->IsendMpiType(&sentNewShadowSimplices[srcIndex][0],
				 sentNewShadowSimplices[srcIndex].size(),
				 mpiType_newShadowSimplices.getType(),
				 ghostExchange.sendRank[srcIndex],
				 &allRequests[4][srcIndex],mpiTagsStart_coarsen+4);
	  }
	
	allRequests[5].resize(ghostExchange.receiveRank.size());
	// send transfered queries (probably mostly empty ...)
#pragma omp parallel for
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    //int index = ghostExchange.sendRank[i];
	    mpiCom->Isend(&queryTransfer[i][0],
			  queryTransfer[i].size(),
			  ghostExchange.receiveRank[i],
			  &allRequests[5][i],mpiTagsStart_coarsen+5);	   
	  }

	DICE_MESH_STATIC std::vector<ShadowSimplex*> transferedReceiveQueries;
	transferedReceiveQueries.clear();
	// now setup the newly received shadow simplices / vertices
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    int count=0;
	    int source=mpiCom->ProbeCount(count,mpiType_newShadowSimplices.getType(),
					  mpiTagsStart_coarsen+4);
	    int srcIndex=receiveGhostRankIndex.at(source);
	    
	    if (count>0) receivedNewShadowSimplices[srcIndex].resize(count);
	    else receivedNewShadowSimplices[srcIndex].clear();	    

	    mpiCom->RecvMpiType(&receivedNewShadowSimplices[srcIndex][0],count,
				mpiType_newShadowSimplices.getType(),source,
				mpiTagsStart_coarsen+4);

	    const std::vector<MpiExchg_NewShadowSimplices> &received=
	      receivedNewShadowSimplices[srcIndex];
	    
	    for (unsigned long j=0;j<received.size();++j)
	      {
		unsigned int index=received[j].getBaseCellIndex();
		GhostSimplex *cur=ghostExchange.receive[source][index];
		//cur->template print<LOG_STD_ALL>();
		
		//if (LocalMesh::ghostSimplexPool.isFree(cur)) exit(0);
		GlobalIdentityValue sGid = received[j].simplexGid;	

		if (sGid == MpiExchg_NewShadowSimplices::neighborIsVoid)
		  {
		    // there is NO shadow neighbor, because the ghost is on a boundary !
		    cur->setNeighbor(cur->cache.c[1],NULL);
		    continue;
		  }
	
		typename GidShadowSimplexHash::iterator sit = 
		  shadowSimplexHash.find(sGid);		
	
		ShadowSimplex *ss;
		
		if (sit==shadowSimplexHash.end())
		  {
		    /*
		    // This is to detect a subtle problem where a merging would result
		    // in the creation of a neighbor shadow simplex although the 
		    // neighbor simplex is actually already a ghost simplex
		    // This problem is due to the fact that shadow simplices do not know
		    // each others ...
		    
		    sit = shadowSimplexHash.find(cur->getGlobalIdentity().get());
		    if (sit != shadowSimplexHash.end())
		      {
			glb::console->print<LOG_ERROR>("FOUND CONFLICT");
			exit(0);
		      }
		    */	

		    int vIndex = received[j].vertexIndex;
		    GlobalIdentityValue vGid = 
		      received[j].vertexGid[vIndex];
		    typename GidShadowVertexHash::iterator vit = 
		      shadowVertexHash.find(vGid);
		    Vertex *v;
		    if (vit==shadowVertexHash.end())
		      {
			ShadowVertex *sv;
			LocalMesh::shadowVertexPool.pop(&sv);
			sv->setCoords(received[j].vertexPos);
			sv->setGeneration(received[j].generation);		
			sv->setData(received[j].vData);			
			sv->setGlobalIdentity(vGid);
			sv->setSetF();
			v=static_cast<Vertex*>(sv);
			shadowVertexHash.insert(std::make_pair(vGid,v));
			// NB: need to do that later, some new shadowVertices may be deleted
			// v->setlocalIndex(shadowVertexHash.getUsedCount()-1);
		      }
		    else v=vit->second;

		    LocalMesh::shadowSimplexPool.pop(&ss);
		    for (int k=0;k<Simplex::NVERT;++k)
		      {
			if (k==vIndex)
			  ss->setVertex(k,v);
			else
			  {
			    Vertex *gv=
			      cur->getVertexByGlobalIdentity(received[j].vertexGid[k]);
			    if (gv==NULL)
			      {
				PRINT_SRC_INFO(LOG_ERROR);
				glb::console->print<LOG_ERROR>
				  ("Error while coarsening, I messed up with the global IDs :( \n");
				cur->template print<LOG_ERROR>();
				ss->template print<LOG_ERROR>();
				glb::console->print<LOG_ERROR>
				  ("When looking for vertex with GID=(%ld,%ld) at baseIndex = %d.\n",
				   (long)GlobalIdentity(received[j].vertexGid[k]).rank(),
				   (long)GlobalIdentity(received[j].vertexGid[k]).id(),
				   received[j].getBaseCellIndex());
				exit(-1);
			      }
			    else ss->setVertex(k,gv);
			  }
		
		      }
		    ss->setNeighbor(vIndex,cur);
		    cur->setNeighbor(cur->cache.c[1],ss);		   
		    ss->setGlobalIdentity(sGid);
		    ss->setGeneration(received[j].simplexGeneration);
		    ss->setData(received[j].sData);	
		    ss->setLocalIndex(LocalMesh::shadowSimplexPool.getUsedCount()-1);
		    ss->setSetF();
		   
		    shadowSimplexHash.insert(std::make_pair(sGid,ss));
		  }
		else
		  {
		    int vIndex = received[j].vertexIndex;
		    ss=static_cast<ShadowSimplex*>(sit->second);		   
		    ss->setNeighbor(vIndex,cur);
		    cur->setNeighbor(cur->cache.c[1],ss);		 
		  }		    	

		// Update shadowExchange.receive only if the query was not transfered
		if (ss->getGlobalIdentity().rank() == source)
		  {		
		    shadowExchange.receive[source].push_back(ss);
		  }
		else
		  {
		    // this is a transfered query, so postpone shadowExchange update !
		    transferedReceiveQueries.push_back(ss);
		  }

	      }
	  }
	
	if (transferedReceiveQueries.size()>0)
	  {
	    // We must order the transfered query so that the remote shadowExchange.send 
	    // and the local shadowExchange.receive arrays match
	    // std::sort(transferedReceiveQueries.begin(),transferedReceiveQueries.end(),
	    //  	      typename ShadowSimplex::cmpPtrLess());

	    // update shadowExchange.receive with transfered queries if necessary
	    for (unsigned long i=0;i<transferedReceiveQueries.size();++i)
	      {
		ShadowSimplex *ss = transferedReceiveQueries[i];
		const int rank=ss->getGlobalIdentity().rank();
		// Check wether we already receive shadowSimplices from that rank and update
		// shadowExchange.receiveRank accordingly !
		// Note that it  may render receiveShadowRankIndex invalid, but we won't 
		// need it anymore
		if (shadowExchange.receive[rank].size()==0)
		  {
		    bool found=0;
		    for (unsigned long j=0;j<shadowExchange.receiveRank.size();++j)
		      if (shadowExchange.receiveRank[j]==rank) found=1;
		    if (!found) 
		      {
			//glb::console->print<LOG_STD_ALL>("ADDING NEW RANK !!!!!!\n");
			shadowExchange.receiveRank.push_back(rank);
			receiveShadowRankIndex.
			  insert(std::make_pair(rank,shadowExchange.receiveRank.size()-1));
		      }
		  }
		shadowExchange.receive[rank].push_back(ss);
	      }
	  }

	std::vector< std::vector<Simplex*> >
	  transferedSendQueries(shadowExchange.sendRank.size());	
	// Receive transfered shadow queries and store them
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  {
	    int count=0;
	    int source=mpiCom->
	      template ProbeCount<GlobalIdentityValue>(count,mpiTagsStart_coarsen+5);
	    
	    DICE_MESH_STATIC std::vector<GlobalIdentityValue> query;	  
	    if (count>0) query.resize(count);
	    else query.clear();
	    //if (count>0) printf("QUERY RECEIVED!\n");
	    //glb::console->print<LOG_STD_ALL>("Received queries size : %d!\n",count);
	    mpiCom->Recv(&query[0],count,source,mpiTagsStart_coarsen+5);
	    if (query.size()>0)
	      {
		const unsigned long nQueries=query.size()/2;
		for (long j=0;j<nQueries;++j)
		  {
		    GlobalIdentity gid(query[2*j]);
		    Simplex *s=ghostExchange.send[source][query[2*j+1]];
		   
		    typename std::map<int,int>::iterator it=
		      sendShadowRankIndex.find(gid.rank());
		    if (it == sendShadowRankIndex.end())
		      {
			//glb::console->print<LOG_STD_ALL>("NEW RANK WAS ADDED !!!!\n");
			shadowExchange.sendRank.push_back(gid.rank());
			transferedSendQueries.resize(transferedSendQueries.size()+1);
			it=sendShadowRankIndex.
			  insert(std::make_pair(gid.rank(),
						shadowExchange.sendRank.size()-1)).
			  first;
		      }
		    transferedSendQueries[it->second].push_back(s);

		    // glb::console->print<LOG_STD_ALL>("ADDING QUERIED SIMPLEX (%ld/%ld): %d->%d->%d\n",j,nQueries,it->first,source,myRank);
		    // s->template print<LOG_STD_ALL>();
		  }
	      }	    
	  }

	// Now we can update shadowExchange.send if necessary
	for (unsigned long i=0;i<shadowExchange.sendRank.size();++i)
	  {
	    
	    if (transferedSendQueries[i].size()>0)
	      {
		// We must order the transfered query so that the remote 
		// shadowExchange.send and the local shadowExchange.send arrays match	
		// std::sort(transferedSendQueries[i].begin(),transferedSendQueries[i].end(),
		// 	  typename Simplex::cmpPtrLocalIndexLess());
		int rank=shadowExchange.sendRank[i];
		shadowExchange.send[rank].insert(shadowExchange.send[rank].end(),
						 transferedSendQueries[i].begin(),
						 transferedSendQueries[i].end());	   
	      }
	  }
	
	unsigned long nmax=0;
	for (long i=0;i<shadowExchange.receiveRank.size();++i)
	  {
	    int rank=shadowExchange.receiveRank[i];
	    if (shadowExchange.receive[rank].size()==0)
	      {
		shadowExchange.reclaimUnusedMemory(shadowExchange.receive[rank]);
		continue;
	      }
		    
	    shadowExchange.receiveRank[nmax++]=rank;

	    std::sort(shadowExchange.receive[rank].begin(),
		      shadowExchange.receive[rank].end(),
		      typename ShadowSimplex::cmpPtrLess());

	    typename std::vector<ShadowSimplex*>::iterator it=
	      std::unique(shadowExchange.receive[rank].begin(),
			  shadowExchange.receive[rank].end());
	    shadowExchange.receive[rank].
	      resize(std::distance(shadowExchange.receive[rank].begin(),it));
	  }
	shadowExchange.receiveRank.resize(nmax);

	nmax=0;
	for (long i=0;i<shadowExchange.sendRank.size();++i)
	  {
	    int rank=shadowExchange.sendRank[i];
	    if (shadowExchange.send[rank].size()==0)
	      {
		shadowExchange.reclaimUnusedMemory(shadowExchange.send[rank]);
		continue;
	      }
	    
	    shadowExchange.sendRank[nmax++]=rank;	
	  
	    std::sort(shadowExchange.send[rank].begin(),
		      shadowExchange.send[rank].end(),
		      typename ShadowSimplex::cmpPtrLocalIndexLess());

	    typename std::vector<Simplex*>::iterator it=
	      std::unique(shadowExchange.send[rank].begin(),
			  shadowExchange.send[rank].end());
	    shadowExchange.send[rank].
	      resize(std::distance(shadowExchange.send[rank].begin(),it));
	  }
	shadowExchange.sendRank.resize(nmax);
     } // nParts>1  
  
    nCoarsened = removedVerticesLID.size();
    long nSimplicesRemoved=removedSimplicesLID.size();
    const unsigned long NSmax=LocalMesh::getNSimplices();
    const unsigned long NVmax=LocalMesh::getNVertices();

    // If coarsening occured, we need to update the local IDs.
    // This is achieved by selecting any simplex/vertex with ID higher
    // than the total number N of simplex/vertex and giving it the index of a 
    // removed simplex/vertex with ID lower than N.
    // FIXME: OPENMP with 2 threads
    if (nCoarsened>0) 
      {
	glb::console->printFlush<LOG_PEDANTIC>("(IDs) ");    

	// First select only those deleted simplices/vertices with
	// ID lower than the current number of simplex/vertex	
	unsigned long count=0;
	for (unsigned long i=0;i<removedSimplicesLID.size();++i)
	  {
	    if (removedSimplicesLID[i]<NSmax)
	      {
		if (i!=count)
		  removedSimplicesLID[count]=removedSimplicesLID[i];
		count++;
	      }
	  }
	removedSimplicesLID.resize(count);
	
	count=0;
	for (unsigned long i=0;i<removedVerticesLID.size();++i)
	  {
	    if (removedVerticesLID[i]<NVmax)
	      {
		if (i!=count)
		  removedVerticesLID[count]=removedVerticesLID[i];
		count++;
	      }
	  }
	removedVerticesLID.resize(count);

	// and update the remaining simplices/vertices local and global IDs where needed
	typename std::vector<LocalIndex>::iterator rms_it = removedSimplicesLID.begin();
	const simplexPtr_iterator its_end=LocalMesh::simplexEnd();
	for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)
	  {
	    if (it->getLocalIndex()>=NSmax) 
	      {
		//glb::console->print<LOG_ERROR>("Changing LOCAL index:");
		//(*it)->template print<LOG_WARNING>();
		it->setLocalIndex(*rms_it); // this also updates the global ID
		//(*it)->template print<LOG_WARNING>();
	
		if ( (++rms_it) == removedSimplicesLID.end()) break;
	      }
	  }

	typename std::vector<LocalIndex>::iterator rmv_it = removedVerticesLID.begin();
	const vertexPtr_iterator itv_end=LocalMesh::vertexEnd();
	for (vertexPtr_iterator it=LocalMesh::vertexBegin();it!=itv_end;++it)
	  {
	    if (it->getLocalIndex()>=NVmax) 
	      {
		it->setLocalIndex(*rmv_it);
		if (it->getGlobalIdentity().rank()==myRank) 
		  {
		    // We may only update globalId if it belongs to the current process
		    it->setGlobalIdentity(myRank,*rmv_it);
		  }
		if ( (++rmv_it) == removedVerticesLID.end()) break;
	      }
	  }
      } // nCoarsened

    // fixme : OPENMP, one task per loop ...
    if (ghostLayerUpdated||(nCoarsenedShared>0)) 
      {		

	// tag all ghost vertices for deletion
	const ghostVertexPtr_iterator itgv_end=LocalMesh::ghostVertexEnd();
	for (ghostVertexPtr_iterator it=LocalMesh::ghostVertexBegin();it!=itgv_end;++it)
	  it->setSetF(false);

	LocalIndex localIndex=0;
	// Free unused ghost simplices, give them a new localID and untag any ghost vertex
	// that belongs to a remaining ghost simplex
	const ghostSimplexPtr_iterator itgs_end=LocalMesh::ghostSimplexEnd();
	for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin();it!=itgs_end;++it)
	  {	    
	    if (!it->isSet()) 
	      {			
		LocalMesh::ghostSimplexPool.recycle(*it);
	      }
	    else
	      {	
		// This is to avoid the (rare) creation of shadow simplices that duplicate
		// ghost simplices in some very particular configurations
		typename GidShadowSimplexHash::iterator sit =
		  shadowSimplexHash.find(it->getGlobalIdentity().get());
		if (sit!=shadowSimplexHash.end())
		  {
		    // When this happens, we must delete the shadow simplex, and correct
		    // the neighborhood by hooking it back to this ghost simplex
		    ShadowSimplex *ss = static_cast<ShadowSimplex*>(sit->second);
		  
		    // tag the shadow simplex for removal
		    ss->setSetF(false);		    
		    // glb::console->print<LOG_WARNING>("REMOVING CONFLICT :\n");
		    // ss->template print<LOG_STD_ALL>();

		    // We tag all the vertices for removal, those that also belong to
		    // the ghost will be untagged anyway
		    for (int i=0;i<Simplex::NVERT;++i)
		      //if (ss->getVertex(i)->isShadow())
			{
			  ss->getVertex(i)->setSetF(false);
			  //ss->getVertex(i)->template print<LOG_STD_ALL>();
			}

		    // reconnect the neighbors
		    for (int i=0;i<Simplex::NNEI;++i)
		      {
			Simplex *nei=ss->getNeighbor(i);
			if (nei!=NULL)
			  {			    
			    int id=nei->getNeighborIndex(ss);
			    nei->setNeighbor(id,*it);
			  }
		      }

		    // So we did reconnect everything but one problem remains now: 
		    // we need to fix the shadowExchange structures, here AND remotely
		  }

		// Now set the new local index and tag the vertices so that they
		// are not removed	
		it->setLocalIndex(localIndex++);
		for (int i=0;i<Simplex::NVERT;++i)
		  it->getVertex(i)->setSetF();
	      }
	  }
	
	localIndex=0;
	// Free unused ghost vertices and give them a new localID 
	const ghostVertexPtr_iterator itgv2_end=LocalMesh::ghostVertexEnd();
	for (ghostVertexPtr_iterator it=LocalMesh::ghostVertexBegin();it!=itgv2_end;++it)
	  if (!it->isSet()) LocalMesh::ghostVertexPool.recycle(*it);
	  else it->setLocalIndex(localIndex++);	
      }
    
    // FIXME: Careful here, this is untested but could be usefull to avoid waisting
    // time sending nothing to processes that are not neighbors anymore !
    if (nParts>1)
      {
	// ghostExchange.template print<LOG_STD_ALL>(" B*ghost simplices*B");
	// shadowExchange.template print<LOG_STD_ALL>("B*shadow simplices*B");
	//bool test=false;
	// Update ghostExchange and free some memory if necessary
	unsigned long newExchgSz=0;
	for (unsigned long i=0;i<ghostExchange.sendRank.size();++i)
	  {
	    int rank=ghostExchange.sendRank[i];
	    if (ghostExchange.send[rank].size()!=0)
	      ghostExchange.sendRank[newExchgSz++]=rank;
	    else
	      ghostExchange.reclaimUnusedMemory(ghostExchange.send[rank]);	   
	  }
	if (newExchgSz!=ghostExchange.sendRank.size())
	  {
	    glb::console->print<LOG_PEDANTIC_ALL>
	      ("INFO: ghostExchange.sendRank has been resized, which means that the number of neighbor processes changed.\n");
	    glb::console->print<LOG_PEDANTIC_ALL>
	      ("This feature is untested, so if I crash, you'll know why ;)");
	    ghostExchange.sendRank.resize(newExchgSz);
	    //test=true;
	  }

	newExchgSz=0;
	for (unsigned long i=0;i<ghostExchange.receiveRank.size();++i)
	  {
	    int rank=ghostExchange.receiveRank[i];
	    if (ghostExchange.receive[rank].size()!=0)
	      ghostExchange.receiveRank[newExchgSz++]=rank;
	    else
	      ghostExchange.reclaimUnusedMemory(ghostExchange.receive[rank]);	   
	  }
	if (newExchgSz!=ghostExchange.receiveRank.size())
	  {
	    glb::console->print<LOG_PEDANTIC_ALL>
	      ("INFO: ghostExchange.receiveRank has been resizes, which means that the number of neighbor processes changed.\n");
	    glb::console->print<LOG_PEDANTIC_ALL>
	      ("This feature is untested, so if I crash, you'll know why ;)");
	    ghostExchange.receiveRank.resize(newExchgSz);
	    //test=true;
	  }
	
	// ghostExchange.template print<LOG_STD_ALL>("ghost simplices");
	// shadowExchange.template print<LOG_STD_ALL>("shadow simplices");
      }
    

    // Deal with the unlikely need to update shadowExchange structure when spurious
    // shadow simplices were created. Most of the time, nothing will be sent !
    // FIXME: OPENMP
    //std::vector<unsigned int> killShadowSend[shadowExchange.receiveRank.size()];
    std::vector< std::vector<unsigned int> > 
      killShadowSend(shadowExchange.receiveRank.size());

    allRequests[6].resize(shadowExchange.receiveRank.size());
    for (unsigned long i=0;i<shadowExchange.receiveRank.size();++i)
      {
	unsigned long newId=0;
	int rank = shadowExchange.receiveRank[i];
	for (unsigned long j=0;j<shadowExchange.receive[rank].size();++j)
	  {
	    if (!shadowExchange.receive[rank][j]->isSet())
	      killShadowSend[i].push_back(j);
	    else 
	      {
		if ( newId != j)
		  shadowExchange.receive[rank][newId]=shadowExchange.receive[rank][j];
		newId++;
	      }
	  }
	if (newId <  shadowExchange.receive[rank].size())
	  shadowExchange.receive[rank].resize(newId);

	mpiCom->Isend(&killShadowSend[i][0],killShadowSend[i].size(),
		      shadowExchange.receiveRank[i],&allRequests[6][i],
		      mpiTagsStart_coarsen+6);
      }


    if (shadowLayerUpdated||(nCoarsenedShared>0)) 
      {
	LocalIndex localIndex=0;
	// Free unused shadow simplices, and check that there are no spurious 
	// shadowVertices remaining. If so, replace them by the corresponding ghost
	const shadowSimplexPtr_iterator itss_end=LocalMesh::shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=LocalMesh::shadowSimplexBegin();it!=itss_end;++it)
	  {
	    if (!it->isSet()) 
	      {
		// This may only happen in the problematic case where a spurious shadow
		// simplex that duplicate an existing ghost was created (see above in 
		// the loop over ghost simplices)
		//exit(-1);
		LocalMesh::shadowSimplexPool.recycle(*it);
	      }
	    else 
	      {
		for (int i=0;i<Simplex::NVERT;++i)
		  {
		    if (!it->getVertex(i)->isSet())
		      {
			GlobalIdentityValue gid=it->getVertex(i)->
			  getGlobalIdentity().get();
			it->setVertex(i,shadowVertexHash.find(gid)->second);
		      }
		  }	
		  
		it->setLocalIndex(localIndex++);		
	      }
	  }	

	localIndex=0;
	// Free unused shadow vertices and update their localId
	const shadowVertexPtr_iterator itsv_end=LocalMesh::shadowVertexEnd();
	for (shadowVertexPtr_iterator it=LocalMesh::shadowVertexBegin();it!=itsv_end;++it)
	  if (!it->isSet()) LocalMesh::shadowVertexPool.recycle(*it);
	  else it->setLocalIndex(localIndex++);
      } // nCoarsened

    // Receive the arrays containing the killed shadowExchange.send simplices
    for (unsigned long i=0;i<shadowExchange.sendRank.size();++i)
      {
	int count=0;
	int source=mpiCom->ProbeCount<unsigned int>(count,mpiTagsStart_coarsen+6);
	
	DICE_MESH_STATIC std::vector<unsigned int> killedShadows;
	if (count>0) killedShadows.resize(count);	
	mpiCom->Recv(&killedShadows[0],count,source,mpiTagsStart_coarsen+6);
	
	if (count>0)
	  {
	    unsigned long newIndex=shadowExchange.send[source].size();
	    for (int j=0;j<count;++j)
	      {
		shadowExchange.send[source][killedShadows[j]]=NULL;
		if (killedShadows[j]<newIndex) newIndex = killedShadows[j];
	      }	    
	    for (unsigned long j=newIndex+1;j<shadowExchange.send[source].size();++j)
	      if (shadowExchange.send[source][j] != NULL)
		shadowExchange.send[source][newIndex++]=
		  shadowExchange.send[source][j];
	    shadowExchange.send[source].resize(newIndex);
	  }
      }
    
     // Now send the new global IDS of the ghostSimplices and their vertices
     std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > 
      sendGlobalIds(ghostExchange.sendRank.size());
    allRequests[7].resize(ghostExchange.sendRank.size());  
#pragma omp parallel for
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];
	sendGlobalIds[i].resize(ghostExchange.send[src].size());
	for (unsigned long j=0;j<sendGlobalIds[i].size();j++)
	  sendGlobalIds[i][j].set(ghostExchange.send[src][j],myRank);

	mpiCom->IsendMpiType(&sendGlobalIds[i][0],sendGlobalIds[i].size(),
			     mpiType_simplicesGlobalIds.getType(),
			     ghostExchange.sendRank[i],&allRequests[7][i],
			     mpiTagsStart_coarsen+7);
      }

    
    // and set the new global Ids ...    
    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int src=mpiCom->Probe(mpiTagsStart_coarsen+7);
	DICE_MESH_STATIC std::vector<MpiExchg_SimplicesGlobalIds> receivedGlobalIds;
	receivedGlobalIds.resize(ghostExchange.receive[src].size());
	mpiCom->RecvMpiType(&receivedGlobalIds[0],receivedGlobalIds.size(),
			    mpiType_simplicesGlobalIds.getType(),src,
			    mpiTagsStart_coarsen+7);
		     
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    GhostSimplex *g=ghostExchange.receive[src][j];
	    MpiExchg_SimplicesGlobalIds &mpiData=receivedGlobalIds[j];
	    for (int k=0;k<Simplex::NVERT;k++)
	      {
		Vertex *v=g->getVertex(k);
		GlobalIdentity globalIdentity(mpiData.vertices[k]);
		// local shared vertices (i.e. those on the interface) 
		// maintain their rank during coarsening, although the ID may change ...
		// NB: this may in some specific case assign a wrong global ID 
		// to a shadow vertex, but this will be corrected when the shadow simplices
		// will be updated.
		if (src == v->getGlobalIdentity().rank())
		  //if (v->getGlobalIdentity().rank()==globalIdentity.rank())
		  v->setGlobalIdentity(globalIdentity);
	      }
	    g->setGlobalIdentity(mpiData.gid);
	  }
      }

    // Same thing for the Shadows
    allRequests[8].resize(shadowExchange.sendRank.size());
    std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > 
      sendShadowGlobalIds(shadowExchange.sendRank.size());
    std::vector<MPI_Request> 
      sentShadowGlobalIdsRequests(shadowExchange.sendRank.size());
#pragma omp parallel for
    for (unsigned long i=0;i<shadowExchange.sendRank.size();i++)
      {
	int src = shadowExchange.sendRank[i];
	sendShadowGlobalIds[i].resize(shadowExchange.send[src].size());
	for (unsigned long j=0;j<sendShadowGlobalIds[i].size();j++)
	  sendShadowGlobalIds[i][j].set(shadowExchange.send[src][j],myRank);

	mpiCom->IsendMpiType(&sendShadowGlobalIds[i][0],sendShadowGlobalIds[i].size(),
			     mpiType_simplicesGlobalIds.getType(),
			     shadowExchange.sendRank[i],&allRequests[8][i],
			     mpiTagsStart_coarsen+8);
      }

    
    // and set the new global Ids ...    
    for (unsigned long i=0;i<shadowExchange.receiveRank.size();i++)
      {
	int src=mpiCom->Probe(mpiTagsStart_coarsen+8);
	DICE_MESH_STATIC std::vector<MpiExchg_SimplicesGlobalIds> receivedGlobalIds;
	receivedGlobalIds.resize(shadowExchange.receive[src].size());
	mpiCom->RecvMpiType(&receivedGlobalIds[0],receivedGlobalIds.size(),
			    mpiType_simplicesGlobalIds.getType(),src,
			    mpiTagsStart_coarsen+8);
		     
	for (unsigned long j=0;j<shadowExchange.receive[src].size();j++)
	  {
	    ShadowSimplex *g=shadowExchange.receive[src][j];
	    MpiExchg_SimplicesGlobalIds &mpiData=receivedGlobalIds[j];
	    for (int k=0;k<Simplex::NVERT;k++)
	      {
		Vertex *v=g->getVertex(k);
		GlobalIdentity globalIdentity(mpiData.vertices[k]);
		// local shared vertices (i.e. those on the interface) 
		// maintain their rank during coarsening, although the ID may change ...
		// NB : checking the source rank here is plain WRONG, don't do it ;)
		//if (src == v->getGlobalIdentity().rank())
		v->setGlobalIdentity(globalIdentity);
	      }
	    g->setGlobalIdentity(mpiData.gid);
	  }
      }

    if (glb::console->willPrint<LOG_DEBUG>())
      {
	shadowExchange.template print<LOG_DEBUG>("SHADOW");
	ghostExchange.template print<LOG_DEBUG>("GHOST");
      }

    // FIXME: remove that (don't need the sum)
    // nCoarsened = mpiCom->sum(nCoarsened);

    if (nCoarsened>0)
      glb::console->print<LOG_INFO>("done. (%ldv/%lds locally removed)\n",
				    nCoarsened,nSimplicesRemoved);
    else
      glb::console->print<LOG_INFO>("done.\n");

    // FIXME: this should be removed ?
    updateCellsCount();    
 
    // Incidence vectors are now invalid !
    if (nCoarsened) incidentSimplices.needFullUpdate=true;

    if ((!glb::console->willPrint<LOG_INFO>())&&(glb::console->willPrint<LOG_STD>())) 
      glb::console->print<LOG_STD>("done.\n");   
   
    if (glb::debug>1)
      dumpToNDnetwork("after_coarsening",
		      IO::NDNET_WithShadows|
		      IO::NDNET_WithGhosts|
		      IO::NDNET_WithNeighbors,
		      true);
    
    // Before returning, we must wait for all the Isends to complete as if we don't
    // their storage might be deallocated before the data is actually sent ...
    for (int i=0;i<nRequestsTotal;++i)
      mpiCom->Waitall(allRequests[i]);

    if (glb::debug) checkConsistencyAndReport<LOG_ERROR>("coarsen");
    
    /*
      {
	bool failed=false;
	//mpiCom->barrier();
	glb::console->printFlush<LOG_STD>("Checking local consistency (coarsen) ... ");  
	if (!LocalMesh::template checkConsistency<LOG_ERROR>(true))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Network is locally inconsistent.\n");
	    failed=true;
	  }
	glb::console->printFlush<LOG_STD>("done.\n"); 
	glb::console->printFlush<LOG_STD>("Checking global consistency (coarsen) ... ");
	if (!checkSharedCellsConsistency<LOG_ERROR>())
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Network is globally inconsistent.\n");
	    failed=true;
	  }
	glb::console->printFlush<LOG_STD>("done.\n"); 

	mpiCom->barrier();
	if (failed) exit(-1);
      }
    */
    return nCoarsened;
  }

  /** \brief Checks the general consistency of the mesh. This is used for debugging purpose.
   */
  template <class LOG>
  bool checkConsistencyAndReport(const std::string &str, bool checkDuplicates=true,
				 bool checkShadows=true, bool checkGlobal=true)
  {
    bool failed=false;
    //mpiCom->barrier();
    if (glb::debug>1)
      {
	dumpToNDnetwork("before_report",
			IO::NDNET_WithShadows|
			IO::NDNET_WithGhosts|
			IO::NDNET_WithNeighbors,
			true);
	mpiCom->barrier();
      }

    glb::console->printFlush<LOG_STD>("Checking local consistency @%s ... ",str.c_str());  
    if (!LocalMesh::template checkConsistency<LOG_ERROR>(checkDuplicates))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Network is locally inconsistent.\n");
	failed=true;
      }
    glb::console->printFlush<LOG_STD>("done.\n"); 
    if (checkGlobal)
      {
	glb::console->printFlush<LOG_STD>("Checking global consistency @%s ... ",
					  str.c_str());
	if (!checkSharedCellsConsistency<LOG_ERROR>(checkShadows))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Network is globally inconsistent.\n");
	    failed=true;
	  }
	glb::console->printFlush<LOG_STD>("done.\n"); 
      }

    mpiCom->barrier();
    if (failed) exit(-1);
    return !failed;
  }

  /** \brief Checks the general consistency of the ghosts and/or shaodw simplices. 
   *  This is used for debugging purpose.
   */
  template <class LOG>
  bool checkSharedCellsConsistency(bool checkShadows=true)
  {
    int myRank=mpiCom->rank();
    bool success=true;

    std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > 
      sendGlobalIds(ghostExchange.sendRank.size());
    std::vector< std::vector<MpiExchg_SimplicesGlobalIds> > 
      receivedGlobalIds;

    /*
    // First check that the rank of any simplex is equal to the rank of its highest
    // ranked vertex.
    // This acutally is only true after a repart or init, this property is not 
    // enforced in genral.
    const simplexPtr_iterator its_end=LocalMesh::simplexEnd();
    for (simplexPtr_iterator it=LocalMesh::simplexBegin();it!=its_end;++it)
      {
	Simplex *s=*it;
	for (int i=0;i<Simplex::NVERT;++i)
	  if (s->getVertex(i)->getGlobalIdentity().rank() > myRank)
	    {
	      glb::console->print<LOG>("Simplex has a vertex with higher rank @%d.\n",i);
	      s->template print<LOG>();
	      s->getVertex(i)->template print<LOG>();
	      success=false;
	    }
      }
    */
    
#pragma omp parallel for
    for (unsigned long i=0;i<ghostExchange.sendRank.size();i++)
      {
	int src = ghostExchange.sendRank[i];
	sendGlobalIds[i].resize(ghostExchange.send[src].size());
	for (unsigned long j=0;j<sendGlobalIds[i].size();j++)
	  sendGlobalIds[i][j].set(ghostExchange.send[src][j],myRank);
      }
   
    ghostExchange.exchangeStruct(sendGlobalIds,receivedGlobalIds,
			   mpiType_simplicesGlobalIds.getType(),true,
			   mpiTagsStart_general+0);

    for (unsigned long i=0;i<ghostExchange.receiveRank.size();i++)
      {
	int fromRank=ghostExchange.receiveRank[i];
	int src = ghostExchange.receiveRank[i];
	for (unsigned long j=0;j<ghostExchange.receive[src].size();j++)
	  {
	    GhostSimplex *g=ghostExchange.receive[src][j];
	    MpiExchg_SimplicesGlobalIds &mpiData=receivedGlobalIds[i][j];
	    if (g->getGlobalIdentity().get() != mpiData.gid)
	      {
		glb::console->print<LOG>("Inconsistent ghost simplex gID found.\n");
		glb::console->print<LOG>("  gID is (%ld,%ld) on remote process %d.\n",
					 (long)GlobalIdentity(mpiData.gid).rank(),
					 (long)GlobalIdentity(mpiData.gid).id(),
					 fromRank);
		glb::console->print<LOG>("  gID is (%ld,%ld) on local process %d.\n",
					 (long)g->getGlobalIdentity().rank(),
					 (long)g->getGlobalIdentity().id(),
					 myRank);
		g->template print<LOG>();		
		success=false;
	      }
	    for (int k=0;k<Simplex::NVERT;k++)
	      {		
		Vertex *v=g->getVertex(k);
		if (mpiData.vertices[k]!=v->getGlobalIdentity().get())
		  {
		    GlobalIdentity gid(mpiData.vertices[k]);
		    glb::console->print<LOG>("Inconsistent vertex gID found.\n");
		    glb::console->print<LOG>("  gID is (%ld,%ld) on remote process %d.\n",
					     (long)gid.rank(),
					     (long)gid.id(),
					     fromRank);
		    glb::console->print<LOG>("  gID is (%ld,%ld) on local process %d.\n",
					     (long)v->getGlobalIdentity().rank(),
					     (long)v->getGlobalIdentity().id(),
					     myRank);	
		    g->template print<LOG>();
		    v->template print<LOG>();
		    success=false;
		  }		
	      }
	    
	  }
      }

    if (!checkShadows) return success;

    sendGlobalIds.clear();
    sendGlobalIds.resize(shadowExchange.sendRank.size());
    receivedGlobalIds.clear();

#pragma omp parallel for
    for (unsigned long i=0;i<shadowExchange.sendRank.size();i++)
      {
	int src = shadowExchange.sendRank[i];
	sendGlobalIds[i].resize(shadowExchange.send[src].size());
	for (unsigned long j=0;j<sendGlobalIds[i].size();j++)
	  sendGlobalIds[i][j].set(shadowExchange.send[src][j],myRank);
      }

    shadowExchange.exchangeStruct(sendGlobalIds,receivedGlobalIds,
			    mpiType_simplicesGlobalIds.getType(),true,
			    mpiTagsStart_general+1);

    for (unsigned long i=0;i<shadowExchange.receiveRank.size();i++)
      {
	int fromRank=shadowExchange.receiveRank[i];
	int src = shadowExchange.receiveRank[i];
	for (unsigned long j=0;j<shadowExchange.receive[src].size();j++)
	  {
	    ShadowSimplex *g=shadowExchange.receive[src][j];
	    MpiExchg_SimplicesGlobalIds &mpiData=receivedGlobalIds[i][j];
	    if (g->getGlobalIdentity().get() != mpiData.gid)
	      {
		glb::console->print<LOG>("Inconsistent shadow simplex gID found.\n");
		glb::console->print<LOG>("  gID is (%ld,%ld) on remote process %d.\n",
					 (long)GlobalIdentity(mpiData.gid).rank(),
					 (long)GlobalIdentity(mpiData.gid).id(),
					 fromRank);
		glb::console->print<LOG>("  gID is (%ld,%ld) on local process %d.\n",
					 (long)g->getGlobalIdentity().rank(),
					 (long)g->getGlobalIdentity().id(),
					 myRank);
		g->template print<LOG>();		
		success=false;
	      }
	    for (int k=0;k<Simplex::NVERT;k++)
	      {		
		Vertex *v=g->getVertex(k);
		if (mpiData.vertices[k]!=v->getGlobalIdentity().get())
		  {
		    GlobalIdentity gid(mpiData.vertices[k]);
		    glb::console->print<LOG>("Inconsistent vertex gID found.\n");
		    glb::console->print<LOG>("  gID is (%ld,%ld) on remote process %d.\n",
					     (long)gid.rank(),
					     (long)gid.id(),
					     fromRank);
		    glb::console->print<LOG>("  gID is (%ld,%ld) on local process %d.\n",
					     (long)v->getGlobalIdentity().rank(),
					     (long)v->getGlobalIdentity().id(),
					     myRank);	
		    g->template print<LOG>();
		    v->template print<LOG>();
		    success=false;
		  }		
	      }
	    
	  }
      }

     return success;
  }

  // virtual const MpiCommunication *getCom() const
  // {
  //   return mpiCom;
  // }

  /** \brief Retrieve the mesh parameters
   * \return the mesh parameters
   */
  const Params &getParams()
  {
    return params;
  }

  /** \brief Retrieve the mesh's bounding box
   *  \param[out] x0 the lower left corner coordinates of the bounding box
   *  \param[out] delta the extent of the bounding box along each dimension
   *  \param world if true, retrieve bounding box dimensions along all NDIM_W dimensions 
   *  (mesh + embedding space. Otherwise, only the NDIM first dimensions are given.
   */
  template <class OutputIterator>
  void getBoundingBox(OutputIterator x0, OutputIterator delta, int world=true)
  {
    if (world)
      {
	std::copy(params.x0,params.x0+NDIM_W,x0);
	std::copy(params.delta,params.delta+NDIM_W,delta);
      }
    else
      {
	std::copy(params.x0,params.x0+NDIM,x0);
	std::copy(params.delta,params.delta+NDIM,delta);
      }
  }
  /*
  // Simple flat simplex gauss integration
  template <class F>
  double quadrature(int nThreads=glb::num_omp_threads, bool global=false)
  {
    typedef GaussQuadratureT<NDIM> GQ;
    double result=0;

#pragma omp parallel num_threads(nThreads) reduction(+:result)
    {
      int th=omp_get_thread_num();
      const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
      for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,nThreads);
	   it!=it_end;++it)
	{
	  F functor(this,*it);
	  result+=GQ::template compute<F>(functor);
	}
    }    
    
    if (global)
      result=mpiCom->sum(result);
    
    return result;
  }
*/

  /** \brief Integrate functor \a f over the mesh
   *  \param nThreads The number of OpenMP threads to use
   *  \param global if false, compute the quadrature for each MPI domain separatly, 
   *  otherwise compute it for the full mesh.
   *  \return the integral of \a f over the mesh.
   */
  template <class F>
  double quadrature(F &f, int nThreads=glb::num_omp_threads, bool global=false)
  {
    typedef GaussQuadratureT<NDIM> GQ;
    double result=0;
   
#pragma omp parallel num_threads(nThreads) reduction(+:result)
    {
      F functor=f;
      int th=omp_get_thread_num();
      const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
      for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,nThreads);
	   it!=it_end;++it)
	{
	  functor.set(this,*it);	  
	  result+=GQ::compute(functor);
	}
    }    
    
    if (global)
      result=mpiCom->sum(result);
    
    return result;
  }

  /** \brief Integrate several functors from a list \a functorList over the mesh
   *  \param result the value of the integral of each functor in the list (as many 
   *  elements as there are functors in the list)
   *  \param nThreads The number of OpenMP threads to use
   *  \param global if false, compute the quadrature for each MPI domai separatly, 
   *  otherwise compute it for the full mesh.   
   */
  template <class FL>
  void quadratureFromList(const FL &functorList, double *result,
			  int nThreads=glb::num_omp_threads, 
			  bool global=false)
  {
    typedef GaussQuadratureT<NDIM> GQ;
    static const int stride = 16; // To avoid false sharing (that's too large, for safety)
    static const int delta = (FL::SIZE+stride);
    std::vector<double> tmpResult(delta*nThreads,0);
    
#pragma omp parallel num_threads(nThreads)
    {
      FL fl = functorList;
      int th=omp_get_thread_num();
      double *res=&tmpResult[delta*th];
      const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
      for (simplexPtr_iterator it=LocalMesh::simplexBegin(th,nThreads);
	   it!=it_end;++it)
	{
	  double tmp[FL::SIZE];
	  internal::mesh::setFunctorList(fl,this,*it);
	  GQ::computeList(fl,tmp);
	  for (int i=0;i<FL::SIZE;++i) res[i]+=tmp[i];
	}
    }    

    for (int i=1;i<nThreads;++i)
      for (int j=0;j<FL::SIZE;++j)
	{
	  tmpResult[j] += tmpResult[i*delta+j];
	}

    std::copy_n(&tmpResult[0],FL::SIZE,result);

    if (global)
      mpiCom->sum(result,FL::SIZE);        
  }

  /** \brief Compute the information necessary to retrieve in constant time the list of
   *  simplices adjacent to any local vertex in the mesh
   *  \param includeGhosts Whether we should include ghost simplices
   *  \param nThreads The number of OpenMP threads to use
   *  \return a pair of a vector of simplices and a vector of indices. Simplices adjacent
   *  to vertex with local index \a V are stored in the range [start,end[ with 
   *  start=&result.first[result.second[V]] and end=&result.first[result.second[V+1]]
   *  \note the two vectors are stored inside the mesh structure so that updates are 
   *  computed only when necessary. Function clearIncidentSimplices() can be used to 
   *  free corresponding memory until next call if required.
   */
  std::pair<const std::vector<Simplex*>&,
	    const std::vector<LocalIndex>&>
  getIncidentSimplices(bool includeGhosts=false,
		       int nThreads=glb::num_omp_threads)
  {
    //#pragma omp critical
    {
      int nNewSimplices=LocalMesh::getNSimplices()-incidentSimplices.oldSimplicesCount;
      bool needFullUpdate=incidentSimplices.needFullUpdate;
     
      // Do a full update if partial update would require too much overhead
      if ( nNewSimplices*10 > LocalMesh::getNSimplices())
	needFullUpdate=true;

      // Forced for now due to a very rare bug in updateIncidentSimplices
      needFullUpdate=true;

      if (needFullUpdate)
	{
	  LocalMesh::getIncidentSimplices
	    (incidentSimplices.simplices,
	     incidentSimplices.index,
	     includeGhosts,nThreads);

	  incidentSimplices.needFullUpdate=false;
	  incidentSimplices.oldSimplicesCount=LocalMesh::getNSimplices();
	  incidentSimplices.oldVerticesCount=LocalMesh::getNVertices();
	  incidentSimplices.oldGhostSimplicesCount=LocalMesh::getNGhostSimplices();
	}
      else if (nNewSimplices>0)
	{
	  LocalMesh::updateIncidentSimplices(incidentSimplices,includeGhosts,nThreads);

	  incidentSimplices.needFullUpdate=false;
	  incidentSimplices.oldSimplicesCount=LocalMesh::getNSimplices();
	  incidentSimplices.oldVerticesCount=LocalMesh::getNVertices();
	  incidentSimplices.oldGhostSimplicesCount=LocalMesh::getNGhostSimplices();
	}      
    }

    return std::pair<const std::vector<Simplex*>&,
		     const std::vector<LocalIndex>&>
      (incidentSimplices.simplices,incidentSimplices.index);    
  }

  /** \brief Free memory allocated getIncidentSimplices to store the lists of simplices
   *  adjacent to each vertex.
   */ 
  void clearIncidentSimplices()
  {
    std::vector<Simplex*> simplices;
    std::vector<LocalIndex> index;   

    // make sure the vectors are deallocated !
    incidentSimplices.simplices.swap(simplices);
    incidentSimplices.index.swap(index);

    incidentSimplices.needFullUpdate=true;
  }
  
protected: 
  typedef typename my_unordered_map<GlobalIdentityValue,void *>::type UMapGlobal;
  typedef typename my_unordered_map<ULong64,void *>::type UMapULL;
  typedef typename my_unordered_map<unsigned long,void *>::type UMapUL;

  typedef typename UMapGlobal::iterator UMapGlobal_iterator;
  typedef typename UMapULL::iterator UMapULL_iterator;
  typedef typename UMapUL::iterator UMapUL_iterator; 

  Params params; // mesh parameters

  MpiCommunication *mpiCom;
  // used for refining
  MpiDataType mpiType_refineQueryResult;
  MpiDataType mpiType_checkRefineResult;   
  MpiDataType mpiType_simplicesGlobalIds;
  MpiDataType mpiType_refineShared;

  // used for coarsening
  MpiDataType mpiType_coarsenShadows;  
  MpiDataType mpiType_coarsenGhosts;
  MpiDataType mpiType_newShadowSimplices;    

  //used for repartitionning
  MpiDataType mpiType_repartSimplices;
  MpiDataType mpiType_repartSimplicesWithCache;
  MpiDataType mpiType_repartVertices;
  MpiDataType mpiType_shadowSimplicesFromGhostQuery;
  MpiDataType mpiType_shadowSimplicesFromGhostReply;    

  int mpiTagsStart_repart;
  int mpiTagsStart_refine;
  int mpiTagsStart_coarsen;
  int mpiTagsStart_general;

  // for shared boundaries communication  
  MpiCellDataExchangeT<Simplex,GhostSimplex> ghostExchange;
  MpiCellDataExchangeT<Simplex,ShadowSimplex> shadowExchange;

  //const Params params;
  //std::vector<double> x0;
  //std::vector<double> delta;  

  //std::vector<unsigned long> globalNCells;
  std::vector<unsigned long> globalNLocalCells[NDIM+1];
  std::vector<unsigned long> globalNCellsCum[NDIM+1];
  double loadImbalanceFactor;

  struct IncidentSimplices {
    typedef LocalIndex IT;

    IncidentSimplices()
    {
      oldSimplicesCount=0;
      oldVerticesCount=0;
      oldGhostSimplicesCount=0;
      needFullUpdate=true;
    }
    
    bool needFullUpdate;
    long oldSimplicesCount;
    long oldVerticesCount;
    long oldGhostSimplicesCount;
    
    std::vector<Simplex*> simplices;
    std::vector<LocalIndex> index;   
  } incidentSimplices;

  template <class Solver, class S>
  void checkRefine_singleSimplex(Solver *solver, S *simplex, std::vector<S *> &cSimplices)
  {
    double result = solver->checkRefine_getValue(simplex);
	      
    if (result<=0) simplex->cache.ptr=NULL; // clean up	      
    else 
      {
	cSimplices.push_back(simplex);
	simplex->cache.pfi.f = result;
#pragma omp task untied
	{
	  double newResult = result;
	  simplex->cache.pfi.i = 
	    solver->checkRefine_getSplitSegmentIndex(simplex,newResult);
	  if (newResult<=0) simplex->cache.ptr=NULL;
	}
      }
    simplex->setTaggedF(false);	
  }

  template <class S>
  void checkRefine_setSimplicesCache_fromCachedSimplices
  (S *solver, std::vector<char> &check, int nThreads,
   std::vector<std::vector<Simplex *> > &cSimplices,
   std::vector<std::vector<GhostSimplex *> > &cGSimplices)
  {
    std::vector<std::vector<Simplex *> > newCSimplices(nThreads);
    std::vector<std::vector<GhostSimplex *> > newCGSimplices(nThreads);

    LocalMesh::resetGhostSimplicesCache(nThreads);
    LocalMesh::resetSimplicesCache(nThreads);
    
#pragma omp parallel num_threads(nThreads)
    {
#pragma omp for nowait
      for (long i=0;i<nThreads;i++)
	{
	  for (auto it=cSimplices[i].begin(); it!=cSimplices[i].end(); ++it)
	    checkRefine_singleSimplex(solver,*it,newCSimplices[i]);
	  for (auto it=cGSimplices[i].begin(); it!=cGSimplices[i].end(); ++it)
	    checkRefine_singleSimplex(solver,*it,newCGSimplices[i]);
	}
    }

    cSimplices.swap(newCSimplices);
    cGSimplices.swap(newCGSimplices);
  }

  template <class S>
  void checkRefine_setSimplicesCache
  (S *solver, std::vector<char> &check, int nThreads,
   std::vector<std::vector<Simplex *> > &cSimplices,
   std::vector<std::vector<GhostSimplex *> > &cGSimplices,
   bool useCache)
  {
    if (useCache) 
      return checkRefine_setSimplicesCache_fromCachedSimplices
	(solver,check,nThreads,cSimplices,cGSimplices);

#pragma omp parallel num_threads(nThreads)
    {
#pragma omp for nowait
      for (long i=0;i<nThreads;i++)
	{
	  cSimplices[i].clear();
	  cGSimplices[i].clear();

	  const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	  for (simplexPtr_iterator it=LocalMesh::simplexBegin(i,nThreads);
	       it!=it_end;++it)
	    checkRefine_singleSimplex(solver,*it,cSimplices[i]);	    
	  
	  if (LocalMesh::getNGhostSimplices()>0)
	    {
	      const ghostSimplexPtr_iterator it_end=LocalMesh::ghostSimplexEnd();
	      for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin(i,nThreads);
		   it!=it_end;++it)
		checkRefine_singleSimplex(solver,*it,cGSimplices[i]);
	    } // if ghostsSimplices

	} //for loop
    } // omp parallel    
  }

  // Eliminates first order conflicts during refinement
  template <class CV>
  void checkRefine_checkConflicts(Simplex *simplex, 
				  std::vector<Simplex*> &cancelled,
				  CV &toRefine)
  {    
    float score = simplex->cache.pfi.f;
    if (score<=0) return;
    int segIndex = simplex->cache.pfi.i;
    // printf("score is %e for simplex %ld (seg %d).\n",
    // 	     score,(long)simplex->getLocalIndex(),segIndex);
	      	      
    SegmentHandle h=simplex->getSegmentHandle(segIndex);
    segment_circulator ci_end=h->getCirculator();
    segment_circulator ci=ci_end;
	      
    // cycle over all simplices containing segment h
    while((++ci) != ci_end)
      {
	bool remove=false;
		    
	//if (fabs(ci->cache.pfi.f-score) <= score*1.E-5)
	if (ci->cache.pfi.f == score)
	  {	
	    // This happens when two simplices want to refine the same simplex
	    // with the exact same score. We define an arbitrary but deterministic
	    // order for them, based on their global identity
	    // NB: the order must be consistent over different processes
		      
	    SegmentHandle otherSeg = ci->getSegmentHandle(ci->cache.pfi.i);
	    if (LocalMesh::compareSegmentHandlesLess(h,otherSeg,true))
	      remove=true;
	    else
	      cancelled.push_back(*ci);
		     
	  }
	else if (ci->cache.pfi.f > score)
	  {
	    // A simplex with higher score already want to refine an edge
	    // of this simplex
	    remove=true;		      
	  }
	else
	  {
	    // A simplex with lower score already want to refine an edge
	    // of this simplex
	    cancelled.push_back(*ci);		      
	  }

	if (remove)
	  {	
	    score=0;
	    break;
	  }
      };
	      
    // now if the score is non-null, we can consider this simplex for
    // refinement.
    if (score>0)
      {
	typedef typename CV::value_type Candidate;
	toRefine.push_back(Candidate(score,h));
      }
    else 
      cancelled.push_back(simplex);   
  }

  // Synchronize boundaries after refining.
  // Used as helper function for refine().
  // sharedVerticesMap maps the globalIdentity of a simplex (before splitting) with the 
  // vertex it created when it was split. sharedVerticesMap only needs to map the newly 
  // created shared vertices and will be updated by the function with the new ghost or 
  // shadow vertices.  
  // The output simplices are sorted such that simplices from 
  // sortedNewSimplices[nNewSimplicesCum[i]] to sortedNewSimplices[nNewSimplicesCum[i+1]-1]
  // will share the new vertex newVertices[i]
  template <class SVM, class VP, class SP, typename TT>
  void synchronizeAfterSplitting(SVM &sharedVerticesMap,
				 VP &newVertexPool, SP &newSimplexPool, 
				 MpiCellDataExchangeT<Simplex,typename SP::Type> 
				 &exchangeInfo,
				 std::vector<Vertex*> &newVertices,
				 std::vector<Simplex*> &sortedNewSimplices,
				 std::vector<TT> &nNewSimplicesCum, int nThreads)
  {   
    typedef typename SVM::iterator SVM_it;
    typedef typename VP::Type SharedVertex;
    typedef typename SP::Type SharedSimplex;

    const int myRank=mpiCom->rank();

    std::vector< std::vector<MpiExchg_RefineShared> > toSend(exchangeInfo.sendRank.size());
    std::vector< std::vector<Simplex*> > newSimplicesToSend(exchangeInfo.sendRank.size());
    
    // setup data for sending
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<exchangeInfo.sendRank.size();i++)
      {	
	const std::vector<Simplex *> &curSend = exchangeInfo.send[exchangeInfo.sendRank[i]];
	if (!MpiExchg_RefineShared::IS_INDEXED) 
	  toSend[i].resize(curSend.size()); // NON INDEXED

	// FIXMEOPENMP not working ...
	//#pragma omp parallel for
	for (unsigned long j=0;j<curSend.size();j++)
	  {
	    Simplex *me=curSend[j];
	 	    
	    if (me->isTagged())
	      {
		// 'me' was split !
		Simplex *other = me->getNeighbor(me->cache.c[0]);
		Vertex *newVertex = me->getVertex(me->cache.c[1]);

		// This is used to rebuild the pointer to the reference simplex, which
		// is partially stored in 'me' and 'other' cache.ui[1]
		union {
		  Simplex *ptr;
		  unsigned int ui[2];
		} tmp;

		tmp.ui[0]=other->cache.ui[1];
		tmp.ui[1]=me->cache.ui[1];
		Simplex *refSimplex = tmp.ptr;
	
		// Simplex *other = me->getPartner();		      
		// Simplex *refSimplex = static_cast<Simplex*>(me->cache.ptr);
		// Vertex *newVertex = static_cast<Vertex*>(other->cache.ptr);
		
		if (MpiExchg_RefineShared::IS_INDEXED)  // INDEXED
		  {
		    MpiExchg_RefineShared res;
		    res.setBaseCellIndex(j);
		    res.set(refSimplex,me,other,newVertex,myRank);
		    toSend[i].push_back(res);
		  }
		else toSend[i][j].set(refSimplex,me,other,newVertex,myRank); // NON INDEXED
		//#pragma omp critical
		newSimplicesToSend[i].push_back(other);
	      }	   
	  }
      }
  
    // we cannot rely on MpiExchg_RefineShared::set() to convert here
    // so we have to deal directly with the MpiStruct
    std::vector< std::vector<typename MpiExchg_RefineShared::MpiStruct> > received;	  
    exchangeInfo.exchangeStruct
      (toSend,received,mpiType_refineShared.getType(),mpiTagsStart_general+2);
    
    // add the new simplices to send in exchangeInfo
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<exchangeInfo.sendRank.size();i++)
      {
	exchangeInfo.send[exchangeInfo.sendRank[i]].
	  insert(exchangeInfo.send[exchangeInfo.sendRank[i]].end(),
		 newSimplicesToSend[i].begin(),
		 newSimplicesToSend[i].end());
      }
 
    typedef typename my_unordered_map<Vertex*,TT>::type VertexPtrMap;
    typedef typename VertexPtrMap::iterator VertexPtrMap_it;
    
    std::vector<Simplex *> newSimplices;
    std::vector<TT> newSimplicesIndex;
    VertexPtrMap vertexMap;
    newVertices.clear();
    
    //FIXME: not sure it is worth parallelizing ...
    for (unsigned long i=0;i<received.size();i++)
      {
	std::vector<SharedSimplex *> &targetSimplexVec = 
	  exchangeInfo.receive[exchangeInfo.receiveRank[i]];
	//#pragma omp parallel for
	for (unsigned long j=0;j<received[i].size();j++)
	  {	    	    
	    typename MpiExchg_RefineShared::MpiStruct *rcvd = &received[i][j];

	    if (rcvd->isEmpty()) continue;

	    Simplex *targetSimplex = 
	      (MpiExchg_RefineShared::IS_INDEXED)?
	      targetSimplexVec[rcvd->getBaseCellIndex()]:
	      targetSimplexVec[j];	   

	    GlobalIdentity newVertexGid(rcvd->newVertex);
	    GlobalIdentity newVertexGeneration(rcvd->newVertexGeneration);
	    Vertex *newVertex;
	    SharedSimplex *newSimplex;	    
	  	   
	    //#pragma omp critical
	    { // start omp critical
	      SVM_it it=sharedVerticesMap.find(rcvd->ref);
	      // Retrieve or create the new vertex :
	      // first check if it was not already created locally : in that case, the key
	      // is the gID of the simplex that caused the split 
	      // -> see LocalMesh::splitSegments
	      if (it==sharedVerticesMap.end())
		{
		  // => create it
		  SharedVertex *sharedVertex;
		  newVertexPool.pop(&sharedVertex);		  
		  sharedVertex->set(rcvd->coords,newVertexPool.getUsedCount()-1);
		  sharedVertex->setData(rcvd->vData);
		  sharedVertex->setGlobalIdentity(newVertexGid);
		  sharedVertex->setGeneration(newVertexGeneration);
		  newVertex=static_cast<Vertex*>(sharedVertex);		  
		  sharedVerticesMap.insert
		    (std::make_pair(rcvd->ref,static_cast<Vertex*>(newVertex)));
		}
	      else
		{		 
		  // => already built, check if its global identity is correct.
		  // former convention was: it must have the rank of the simplex that
		  // caused the split, but this does not work for some ghost vertices
		  // Convention: the gID with lowest rank wins ...
		  newVertex=static_cast<Vertex*>(it->second);	
		  // if (newVertex->getGlobalIdentity().rank() != 
		  //     GlobalIdentity(rcvd->ref).rank())
		  if (newVertex->getGlobalIdentity().rank() >
		      newVertexGid.rank())
		    {	
		      newVertex->setGlobalIdentity(newVertexGid);	
		      newVertex->setGeneration(newVertexGeneration);	
		    }
		}

	      // create the new simplex
	      newSimplexPool.pop(&newSimplex);
	      newSimplex->setLocalIndex(newSimplexPool.getUsedCount()-1);
	      newSimplex->setGlobalIdentity(rcvd->other);
	      newSimplex->setGeneration(rcvd->otherGeneration);
	   
	    } // end omp critical

	    Vertex *v0=targetSimplex->getVertexByGlobalIdentity(rcvd->segment[0]);
	    Vertex *v1=targetSimplex->getVertexByGlobalIdentity(rcvd->segment[1]);

	    if ((v0==NULL)||(v1==NULL))
	      {
		targetSimplex->template print<LOG_ERROR>("\n");
		if (v0!=NULL)
		  v0->template print<LOG_ERROR>();
		else 
		  {
		    GlobalIdentity gid(rcvd->segment[0]);
		    glb::console->
		      template print<LOG_ERROR>("v0 is NULL (gid=(%ld,%ld))\n",
						(long)gid.rank(),(long)gid.id());
		  }
		
		if (v1!=NULL)
		  v1->template print<LOG_ERROR>();
		else 
		  {
		    GlobalIdentity gid(rcvd->segment[1]);
		    glb::console->
		      template print<LOG_ERROR>("v1 is NULL (gid=(%ld,%ld))\n",
						(long)gid.rank(),(long)gid.id());
		  }
		PRINT_SRC_INFO(LOG_ERROR);
		glb::console->print<LOG_ERROR>
		  ("Vertices global identities seem to be out of sync !\n");
		exit(-1);
	      }

	    SegmentHandle segment(Segment(v0,v1,NULL));
	    targetSimplex->splitSegment(segment,NULL,newVertex,newSimplex);
	    
	    targetSimplex->setData(rcvd->sData[0]);
	    newSimplex->setData(rcvd->sData[1]);

	    //#pragma omp critical (SYNC_PB)
	    {
	      targetSimplexVec.push_back(newSimplex);	 
	      newSimplices.push_back(static_cast<Simplex*>(targetSimplex));
	      std::pair<VertexPtrMap_it,bool> tmp=
		vertexMap.insert(std::make_pair(newVertex,vertexMap.size()));
	      if (tmp.second) newVertices.push_back(newVertex);
	      TT index = tmp.first->second;
	      newSimplicesIndex.push_back(index);
	    }
	  }
      }
    //if (global_stop_signal) {mpiCom->barrier();exit(0);}

    // prepare vectors of newly created simplices and vertices
    // to be used by the caller to fix neighbors ...
    nNewSimplicesCum.assign(vertexMap.size()+1,0);
    sortedNewSimplices.resize(newSimplices.size());
    for (unsigned long i=0;i<newSimplicesIndex.size();i++)
      nNewSimplicesCum[newSimplicesIndex[i]+1]++;
    for (unsigned long i=1;i<nNewSimplicesCum.size();i++)
      nNewSimplicesCum[i]+=nNewSimplicesCum[i-1];
    std::vector<TT> nNewSimplicesCumTmp=nNewSimplicesCum;
    for (unsigned long i=0;i<newSimplicesIndex.size();i++)
      sortedNewSimplices[nNewSimplicesCumTmp[newSimplicesIndex[i]]++]=newSimplices[i];
    
    // and update the cumulated counts in exchangeInfo    
    exchangeInfo.updateNCum();
  }
  
     
  // Builds a tesselation from a regular grid and a partition of its NDIM-simplices.
  // The partition contains the indices of the NDIM-simplices that belong to the local node
  // or may be empty (i.e. all nodes have all the data)
  //template <class SG, typename VT>
  //void initFromRegularGrid(SG *sg, const std::vector<VT> &partition)
  template <class IST, typename VT>
  void initFromImplicitTesselation(IST *st, const std::vector<VT> &partition)
  {
    typedef typename IST::Cell ITCell;
    typedef typename my_unordered_map<VT,void*>::type UMapVT; 
    typedef typename UMapVT::iterator UMapVT_iterator;

    const bool usePartition = (partition.size()==0)?false:true;  
    const int nParts = mpiCom->size();
    const int myRank = mpiCom->rank();

    // number of cells over all tesselation
    std::vector<unsigned long> nCellsTotal=st->getNCells();
    // number of cells per voxel (i.e. grid volume element)
    std::vector<unsigned long> nCellsPerVox=st->getNCellsPerVoxel(); 
    // local number of cells
    std::vector<unsigned long> nCells(NDIM+1); 
        
    if (sizeof(ULong64)!=8)
      {
	if (sizeof(ULong64)<8)
	  {
	    PRINT_SRC_INFO(LOG_WARNING);
	    glb::console->print<LOG_WARNING>
	      ("Type ULong64 is smaller than expected. (%d bytes, should be %d)\n",
	       sizeof(ULong64),8);
	    glb::console->print<LOG_WARNING>("You should edit file 'types.hxx'.\n");
	  }
	else
	  {
	    glb::console->print<LOG_STD>
	      ("WARNING: Type ULong64 is bigger than expected. (%d bytes, should be %d)\n",
	       sizeof(T),8);
	    glb::console->print<LOG_STD>("  You should edit file 'types.hxx'.\n");
	  }
      }
    
    if (!usePartition)
      {
	nCells=nCellsTotal;
      }
    else
      {
	nCells[NDIM]=partition.size();
	// here, nCells[0...NDIM-1] are estimates only
	unsigned long nvox=nCells[NDIM]/nCellsPerVox[NDIM] + 1; 
	for (int i=0;i<NDIM;i++)
	  {
	    nCells[i]=nCellsPerVox[i]*nvox;
	  }	
	
	glb::console->printToBuffer<LOG_PEDANTIC>
	  ("Estimated local cells count (node 0): [%ld",nCells[0]);
	for (int i=1;i<=NDIM;i++) 
	  glb::console->printToBuffer<LOG_PEDANTIC>(",%ld",nCells[i]);
	glb::console->printToBuffer<LOG_PEDANTIC>("] k-cells.\n");
	glb::console->flushBuffer<LOG_PEDANTIC>();   	
      }  
   
    glb::console->print<LOG_PEDANTIC>
      ("Expected structure size per volume element: %ld x %ld + %ld x %ld = %ld Bytes.\n",
       nCellsPerVox[0],sizeof(Vertex),nCellsPerVox[NDIM],sizeof(Simplex),  
       nCellsPerVox[0]*sizeof(Vertex)+nCellsPerVox[NDIM]*sizeof(Simplex));  
					  		

    // estimate of the number of shadow and ghost simplices for the pools
    long nShadowEstimate=0; 
    long nGhostEstimate=0;
    if (usePartition)
      {
	nShadowEstimate=pow(nCellsTotal[NDIM]/nParts,((double)NDIM-1)/((double)NDIM))*3.0;
	nGhostEstimate=pow(nCellsTotal[NDIM]/nParts,((double)NDIM-1)/((double)NDIM))*3.0;
	glb::console->print<LOG_PEDANTIC>
	  ("Estimated ghost/shadow simplices count : %ld/%ld.\n",
	   nGhostEstimate,nShadowEstimate);
      }
    //nGhostEstimate=0;

    glb::console->printFlush<LOG_STD>("Creating simplices ... ");
   
    LocalMesh::vertexPool.reserve(nCells[0]);    
    LocalMesh::simplexPool.reserve(nCells[NDIM]); 
    if (nShadowEstimate) 
      LocalMesh::shadowSimplexPool.reserve(nShadowEstimate);
    if (nGhostEstimate) 
      LocalMesh::ghostSimplexPool.reserve(nGhostEstimate);
    /*
    LocalMesh::vertexPool.setAllocChunkSize(params.allocFactor*nCells[0]);
    LocalMesh::simplexPool.setAllocChunkSize(params.allocFactor*nCells[NDIM]);
    if (nShadowEstimate) 
      LocalMesh::shadowSimplexPool.setAllocChunkSize(params.allocFactor*nShadowEstimate);
    if (nGhostEstimate) 
      LocalMesh::ghostSimplexPool.setAllocChunkSize(params.allocFactor*nGhostEstimate);
    */
    // Vertices storage
    std::vector<Vertex *> vertexArr; // used when there is no partition, faster ...
    UMapUL vertexTable; // used only when partitioning
    UMapUL shadowVertexTable; // used only when partitioning
    UMapUL ghostVertexTable; // used only when partitioning
    // Simplices storage
    std::vector<Simplex *> simplexArr;
    UMapVT simplexTable;    
    UMapVT shadowSimplexTable;  
    UMapVT ghostSimplexTable; 
    if (usePartition)
      {
	// reserve is not available on several implementations ...
	vertexTable.rehash(nCells[0]/
			   (vertexTable.max_load_factor()-0.01)); 	
	shadowVertexTable.rehash(nShadowEstimate/
				 (shadowVertexTable.max_load_factor()-0.01));
	ghostVertexTable.rehash(nGhostEstimate/
				(ghostVertexTable.max_load_factor()-0.01));

	simplexTable.rehash(nCells[NDIM]);
	shadowSimplexTable.rehash(nShadowEstimate/
				  (shadowSimplexTable.max_load_factor()-0.01));
	ghostSimplexTable.rehash(nGhostEstimate/
				 (ghostSimplexTable.max_load_factor()-0.01));
      }
    else
      {
	vertexArr.resize(nCells[0]);  
	for (unsigned long i=0;i<nCells[0];i++) 
	  {
	    Coord pos[NDIM_W]={0};
	    //std::vector<Coord> pos(NDIM_W,0);
	    LocalMesh::vertexPool.pop(&vertexArr[i]);
	    st->getPosition(ITCell(0,i),&pos[0]);
	    vertexArr[i]->set(&pos[0],i);
	    vertexArr[i]->setGeneration(0,i);
	    vertexArr[i]->setGlobalIdentity(mpiCom->rank(),i);
	  }    

	simplexArr.resize(nCells[NDIM]);
	for (unsigned long i=0;i<nCells[NDIM];i++) 
	  {
	    LocalMesh::simplexPool.pop(&simplexArr[i]);
	    simplexArr[i]->setLocalIndex(i);
	  }
      }

    // Get a projected volume functor
    const typename LocalMesh::SimplexFunctor *projectedVolumePtr = 
      LocalMesh::getSimplexFunctorPtr("signedProjectedVolume");

    if (projectedVolumePtr == NULL)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("Could not find 'projectedVolume' functor for simplices.\n");	
	exit(-1);
      }
    const typename LocalMesh::SimplexFunctor &projectedVolumeFunctor=(*projectedVolumePtr);
       
    // Create each simplex, one by one ...
    for (unsigned long i=0;i<nCells[NDIM];i++)
      {
	unsigned long index=(usePartition)?partition[i]:i;
	Vertex *v[Simplex::NVERT];
	//std::vector<ITCell> cellV;
	ITCell cellV[Simplex::NVERT];
	
	// create or retrieve the NDIM+1 vertices
	//st->getVertices(ITCell(NDIM,index),std::back_inserter(cellV));
	st->getVertices(ITCell(NDIM,index),&cellV[0]);
	
	if (usePartition)
	  {
	    for (long j=0;j<Simplex::NVERT;j++) 
	      {
		unsigned long id = cellV[j].id();
		UMapUL_iterator it = vertexTable.find(id);
		if (it==vertexTable.end())
		  {
		    Coord pos[NDIM_W]={0};
		    //std::vector<Coord> pos(NDIM_W,0);
		    LocalMesh::vertexPool.pop(&v[j]);
		    st->getPosition(ITCell(0,id),&pos[0]);
		    v[j]->set(&pos[0],LocalMesh::vertexPool.getUsedCount()-1);
		    v[j]->setGlobalIdentity(mpiCom->rank(),v[j]->getLocalIndex());
		    v[j]->setGeneration(0,id);
		    vertexTable.insert
		      (std::make_pair(id,(void*)static_cast<Vertex*>(v[j])));
		  }
		else v[j]=static_cast<Vertex *>(it->second);		
	      }
	    
	    // set the configuration of the current simplex
	    // the simplex neighbors do not have to be set yet
	    Simplex *simplex;
	    LocalMesh::simplexPool.pop(&simplex);
	    simplex->setElements(v,projectedVolumeFunctor);
	    simplex->setLocalIndex(LocalMesh::simplexPool.getUsedCount()-1);
	    simplex->setGeneration(0,index);
	    simplexTable.insert
	      (std::make_pair((VT)index,(void*)static_cast<Simplex*>(simplex)));
	  }
	else
	  {
	    for (long j=0;j<Simplex::NVERT;j++) 
	      {
		v[j]=vertexArr[cellV[j].id()];
		v[j]->setGeneration(0,cellV[j].id());
	      }

	    // set the configuration of the current simplex
	    // the simplex neighbors do not have to be set yet
	    Simplex *simplex=simplexArr[i];
	    simplex->setElements(v,projectedVolumeFunctor);
	    simplex->setGeneration(0,index);
	  }
      } 

    glb::console->printFlush<LOG_STD>("done.\n");
  
    // If needed, build the ghost layer
    if (nGhostEstimate)
      {
	glb::console->printFlush<LOG_STD>("Building ghost layer ... ");
	std::set<VT> boundaryVertex;
	
	// first retrieve all the boundary vertices
	glb::console->printFlush<LOG_PEDANTIC>("(boundary) ");
#pragma omp parallel for
	for (unsigned long i=0;i<nCells[NDIM];i++)
	  {
	    std::vector<ITCell> cellV(Simplex::NNEI,ITCell::empty);
	    st->getNeighbors(ITCell(NDIM,partition[i]),&cellV[0]);

	    for (int j=0;j<Simplex::NNEI;j++) 
	      {
		if (cellV[j].isEmpty()) continue;
		VT neiId=cellV[j].id();
		UMapVT_iterator it = simplexTable.find(neiId);
		if (it==simplexTable.end())
		  {
		    //std::vector<ITCell> vertex(Simplex::NVERT);
		    ITCell vertex[Simplex::NVERT];

		    st->getVertices(ITCell(NDIM,neiId),&vertex[0]);
		    for (int k=0;k<Simplex::NVERT;k++) 
		      {
			UMapUL_iterator it = vertexTable.find(vertex[k].id());
			if (it!=vertexTable.end())
			  {
#pragma omp critical
			    boundaryVertex.insert((VT)vertex[k].id());
			  }
		      }
		  }
	      }
	  }
	
	// ghost simplices are non local simplices with at least one vertex in the  
	// local distribution
	glb::console->printFlush<LOG_PEDANTIC>("(create) ");
	for (typename std::set<VT>::iterator it=boundaryVertex.begin();
	     it!=boundaryVertex.end();++it)
	  {
	    std::vector<ITCell> cellV;		
	    st->getNeighbors(ITCell(0,*it),cellV,NDIM);
	
	    for (unsigned long i=0;i<cellV.size();i++)
	      {
		VT neiId = cellV[i].id();	
		UMapVT_iterator sit = simplexTable.find(neiId);

		if (sit==simplexTable.end())
		  {
		    sit = ghostSimplexTable.find(neiId);
		    if (sit==ghostSimplexTable.end())
		      {
			GhostSimplex *ghost;
			LocalMesh::ghostSimplexPool.pop(&ghost);
			//std::vector<ITCell> vertex(Simplex::NVERT);
			ITCell vertex[Simplex::NVERT];
			Vertex *v[Simplex::NVERT];
			st->getVertices(ITCell(NDIM,neiId),&vertex[0]);

			for (int j=0;j<Simplex::NVERT;j++)
			  {
			    UMapUL_iterator vit = vertexTable.find(vertex[j].id());
			    if (vit==vertexTable.end())
			      {
				vit = ghostVertexTable.find(vertex[j].id());
				if (vit==ghostVertexTable.end())
				  {	
				    GhostVertex *tmpV;
				    //std::vector<Coord> pos(NDIM_W,0);
				    Coord pos[NDIM_W]={0};
				    LocalMesh::ghostVertexPool.pop(&tmpV);
				    st->getPosition(vertex[j],&pos[0]);
				    tmpV->set(&pos[0],LocalMesh::ghostVertexPool.getUsedCount()-1);
				    ghostVertexTable.
				      insert
				      (std::make_pair
				       ((unsigned long)vertex[j].id(),
					(void*)static_cast<Vertex*>(tmpV)));

				    v[j]=tmpV;
				  } else v[j]=static_cast<Vertex*>(vit->second);
			      } else v[j]=static_cast<Vertex*>(vit->second);
			    
			    v[j]->setSharedF(true);
			    v[j]->setGeneration(0,vertex[j].id());
			    v[j]->setGlobalIdentity(GlobalIdentity::EMPTY_RANK,
						    vertex[j].id());
			  }
			
			ghost->setElements(v,projectedVolumeFunctor);
			ghost->setLocalIndex
			  (LocalMesh::ghostSimplexPool.getUsedCount()-1);

			ghost->setGlobalIdentity(GlobalIdentity::EMPTY_RANK,neiId);
			ghost->setGeneration(0,neiId);
			ghostSimplexTable.insert
			  (std::make_pair(neiId,
					  (void*)static_cast<Simplex*>(ghost)));
		      }
		  }
	      }
	  }

	if (glb::console->willPrint<LOG_INFO>())
	  glb::console->print<LOG_INFO>
	    ("done. (+%ld/%ld ghost simplices/vertices)\n",
	     LocalMesh::ghostSimplexPool.getUsedCount(),
	     LocalMesh::ghostVertexPool.getUsedCount());

	else glb::console->print<LOG_STD>("done.\n");
      }

    mpiCom->barrier();
    //exit(0);
    glb::console->printFlush<LOG_STD>("Identifying neighbors ... ");

    // Set the simplices neighbors
#pragma omp parallel for
    for (unsigned long i=0;i<nCells[NDIM];i++)
      {
	unsigned long index=(usePartition)?partition[i]:i;	
	Simplex *nei[Simplex::NNEI];
	std::vector<ITCell> cellV(Simplex::NNEI,ITCell::empty);
	Simplex *simplex;

	// retrieve the NDIM+1 neighboring simplices
	st->getNeighbors(ITCell(NDIM,index),&cellV[0]);

	if (usePartition)
	  {
	    simplex=static_cast<Simplex*>(simplexTable.find(index)->second);
	    for (long j=0;j<Simplex::NVERT;j++) 
	      {	
		if (cellV[j].isEmpty()) 
		  {
		    // simplex touches the boundary
		    nei[j]=static_cast<Simplex*>(NULL);
		    continue;
		  }	     	
		  
		// Check if the neighbor exists
		VT neiId=cellV[j].id();
		bool found=false;
		UMapVT_iterator it = simplexTable.find(neiId);
		if (it==simplexTable.end())
		  {
		    it = ghostSimplexTable.find(neiId);		
		    if (it != ghostSimplexTable.end()) found=true;
		  }
		else found=true;
		  
		// if not, it is a shadow simplex
		if (!found)
		  {
		    std::pair<ShadowSimplex *, Vertex *> result; 
#pragma omp critical
		    result = findOrCreateShadowSimplex(st,projectedVolumeFunctor,
						       neiId,index,
						       shadowSimplexTable,vertexTable,
						       ghostVertexTable,shadowVertexTable);
		    ShadowSimplex *shadow = result.first;
		    Vertex *oppVertex = result.second;		   
		 
		    shadow->addNeighbor(oppVertex,simplex);	
		    nei[j]=shadow;
		  }
		else nei[j]=static_cast<Simplex*>(it->second); // regular local neighbor, nothing special		
	      }
	  }
	else // no partition => no shadow
	  {	    
	    simplex=simplexArr[i];
	    for (long j=0;j<Simplex::NNEI;j++) 
	      {
		if (cellV[j].isEmpty()) // simplex touches the boundary
		  nei[j]=NULL;
		else
		  nei[j]=simplexArr[cellV[j].id()];		    
	      }
	  }
	
	//#pragma omp critical
	simplex->setNeighbors(nei);
      }
    
    // and, if necessary, set the ghost simplices neighbors
    if (nGhostEstimate)
      {
	glb::console->printFlush<LOG_PEDANTIC>("(ghosts) ");

#pragma omp parallel for
	for (long i=0;i<glb::num_omp_threads;i++)
	  {	
	    const ghostSimplexPtr_iterator it_end=LocalMesh::ghostSimplexEnd();
	    for (ghostSimplexPtr_iterator it=LocalMesh::ghostSimplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	      {
		GhostSimplex *simplex=*it;
		Simplex *nei[Simplex::NNEI];
		unsigned long index=simplex->getGlobalIdentity().id();
		std::vector<ITCell> cellV(Simplex::NNEI,ITCell::empty);
		st->getNeighbors(ITCell(NDIM,index),&cellV[0]);
		for (long j=0;j<Simplex::NVERT;j++) 
		  {	
		    if (cellV[j].isEmpty()) 
		      {
			// simplex touches the boundary
			nei[j]=(Simplex*)NULL;
			continue;
		      }	
		    // Check if the neighbor exists
		    VT neiId=cellV[j].id();
		    bool found=false;
		    UMapVT_iterator it = simplexTable.find(neiId);
		    if (it==simplexTable.end())
		      {
			it = ghostSimplexTable.find(neiId);		
			if (it != ghostSimplexTable.end()) found=true;
		      }
		    else found=true;

		    // if not, it is a shadow simplex
		    if (!found)
		      {
			std::pair<ShadowSimplex *, Vertex *> result; 
#pragma omp critical
			result = findOrCreateShadowSimplex(st,projectedVolumeFunctor,
							   neiId,index,
							   shadowSimplexTable,vertexTable,
							   ghostVertexTable,
							   shadowVertexTable);
			ShadowSimplex *shadow = result.first;
			Vertex *oppVertex = result.second;
			shadow->addNeighbor(oppVertex,simplex);	
			//shadow->addNeighbor((Vertex*)itV->second,simplex);	
			nei[j]=shadow;
		      }
		    else nei[j]=static_cast<Simplex*>(it->second); // regular local neighbor, nothing special		
		  }
#pragma omp critical
		simplex->setNeighbors(nei);	
	      }	   
	  }
      }
    
    // set correct allocation factor for structure ...
    //LocalMesh::updatePoolsAllocSize();
    /*
    nCells=LocalMesh::getNCells();
    LocalMesh::vertexPool.setAllocChunkSize(params.allocFactor*LocalMesh::getNVertices());
    LocalMesh::simplexPool.setAllocChunkSize(params.allocFactor*LocalMesh::getNSimplices());
    if (nShadowEstimate) 
      LocalMesh::shadowSimplexPool.setAllocChunkSize(params.allocFactor*LocalMesh::getNShadowSimplices()); 
    if (nGhostEstimate) 
      LocalMesh::ghostSimplexPool.setAllocChunkSize(params.allocFactor*LocalMesh::getNGhostSimplices()); 
    */
    glb::console->print<LOG_STD>("done.\n");
    
    if (usePartition)
      {
	// print some info
	if (glb::console->willPrint<LOG_PEDANTIC>())
	  {	    
	    std::vector<unsigned long> nCellsFinalMax=LocalMesh::getNCells();
	    std::vector<unsigned long> nCellsFinalMin=LocalMesh::getNCells();
	    mpiCom->max(&nCellsFinalMax[0],nCellsFinalMax.size());
	    mpiCom->min(&nCellsFinalMin[0],nCellsFinalMin.size());
	    glb::console->printToBuffer<LOG_PEDANTIC>("Actual local cells count (min/max): [%ld/%ld",
						 nCellsFinalMin[0],nCellsFinalMax[0]);
	    for (int i=1;i<=NDIM;i++) glb::console->printToBuffer<LOG_PEDANTIC>(",%ld/%ld",
									   nCellsFinalMin[i],nCellsFinalMax[i]);
	    glb::console->printToBuffer<LOG_PEDANTIC>("] k-cells.\n");
	    glb::console->flushBuffer<LOG_PEDANTIC>();  
	  }
	
	unsigned long nShadowMax=mpiCom->max(LocalMesh::getNShadowSimplices());
	unsigned long nShadowMin=mpiCom->min(LocalMesh::getNShadowSimplices());
	glb::console->print<LOG_INFO>("Number of shadow simplices: %ld / %ld (min/max).\n",nShadowMin,nShadowMax);

	unsigned long nShadowVMax=mpiCom->max(LocalMesh::getNShadowVertices());
	unsigned long nShadowVMin=mpiCom->min(LocalMesh::getNShadowVertices());
	glb::console->print<LOG_INFO>("Number of shadow vertices: %ld / %ld (min/max).\n",nShadowVMin,nShadowVMax);
     	
	glb::console->printFlush<LOG_STD>("Setting up shadow cells ... ");
	
	// retrieve vectors of pointers to all shadow/ghost simplices
	glb::console->printFlush<LOG_PEDANTIC>("(retrieve) ");	
	std::vector<GhostSimplex *> ghostArr = LocalMesh::getGhostSimplicesArray();
	std::vector<ShadowSimplex *> shadowArr = LocalMesh::getShadowSimplicesArray();

	// create the tables that map local and remote shadow/ghost simplices
	glb::console->printFlush<LOG_PEDANTIC>("(com) ");
	initExchangeCells(ghostArr,simplexTable,vertexTable,ghostExchange);	
	initExchangeCells(shadowArr,simplexTable,vertexTable,shadowExchange);
	//shadowSimplicesSend = shadowExchange.send;
	//shadowSimplicesReceive = shadowExchange.receive;
	//sendNeighborNodesRank = shadowExchange.sendRank;
	//receiveNeighborNodesRank = shadowExchange.receiveRank;
	glb::console->print<LOG_STD>("done.\n");	  

	//printExchangeData<LOG_DEBUG>(ghostExchange,"ghost cells");
	//printExchangeData<LOG_DEBUG>(shadowExchange,"shadow cells");
	ghostExchange.template print<LOG_DEBUG>("ghost cells");
	shadowExchange.template print<LOG_DEBUG>("shadow cells");

	// get a correct globalIdentity for all shared vertices
	// FIXME: use alltoall communication !
	glb::console->printFlush<LOG_STD>("Setting up shared vertices ... ");
	std::vector<Vertex*> sharedVerticesArr = LocalMesh::getSharedVerticesArray();
	for (int i=0;i<nParts;i++)
	  {
	    std::vector<GlobalIdentityValue> identityArr;
	    unsigned long len=sharedVerticesArr.size();

	    if (myRank==i)
	      {	
		identityArr.resize(len);
		for (unsigned long j=0;j<len;j++)
		  identityArr[j] = sharedVerticesArr[j]->getGlobalIdentity().id();	
	      }
	    
	    mpiCom->Bcast(&len,i,1);
	    
	    if (myRank!=i) identityArr.resize(len);	    

	    mpiCom->Bcast(identityArr,i);
	    
	    for (unsigned long j=0;j<len;j++)
	      {
		UMapUL_iterator it = vertexTable.find(identityArr[j]);
		if (it==vertexTable.end())
		  {
		    identityArr[j]=GlobalIdentity::max.get();
		    //identityArr[j]=GlobalIdentity(0,0).get();
		  }
		else
		  {
		    Vertex *v=static_cast<Vertex*>(it->second);
		    identityArr[j]=GlobalIdentity(myRank,v->getLocalIndex()).get();
		  }
	      }
	    // FIXME: this is not optimal, nodes with lower rank may get more important charge ...
	    mpiCom->Allreduce_inplace(identityArr,MPI_MIN); 

	    if (myRank==i)
	      {
		for (unsigned long j=0;j<len;j++)
		  sharedVerticesArr[j]->setGlobalIdentity(identityArr[j]);
	      }	  
	  }
	glb::console->printFlush<LOG_STD>("done.\n");
      }

    glb::console->printFlush<LOG_STD>("Initializing tree structure ... ");
    Tree::initRoot(LocalMesh::simplexBegin(),LocalMesh::simplexEnd(),LocalMesh::getNCells(NDIM));
    //mpiCom->barrier();exit(0);
    updateCellsCount();
    glb::console->print<LOG_STD>("done.\n");  
    LocalMesh::initialized=true; 
  }
  
  // This function retrieves or creates a shadow simplex from the simplex with id 
  // shadowId and neighbor neigborId
  // This is used as a helper for initFromRegularGrid(...)
  template <class IST,class ST, class VT, typename IDTA, typename IDTB, class PVF>
  std::pair<ShadowSimplex *, Vertex *> findOrCreateShadowSimplex
  (IST *st,
   const PVF &projectedVolumeFunctor,
   IDTA shadowId, 
   IDTB neighborId,
   ST &shadowSimplexTable, 
   VT &vertexTable,
   VT &ghostVertexTable, 
   VT &shadowVertexTable)
  {  
    typedef typename IST::Cell ITCell;   
    ShadowSimplex *shadow;
    std::vector<ITCell> shadowVert(Simplex::NVERT);
    std::vector<ITCell> neighborVert(Simplex::NVERT);
    typename ST::iterator itS=shadowSimplexTable.find(shadowId);

    st->getVertices(ITCell(NDIM,shadowId),&shadowVert[0]);
    st->getVertices(ITCell(NDIM,neighborId),&neighborVert[0]);

    if (itS==shadowSimplexTable.end())
      {
	LocalMesh::shadowSimplexPool.pop(&shadow);
	//ShadowVertex *v[Simplex::NVERT];
	Vertex *v[Simplex::NVERT];
	for (int j=0;j<Simplex::NVERT;j++)
	  {
	    unsigned long id = shadowVert[j].id();
	    typename VT::iterator it = vertexTable.find(id);
	    if (it==vertexTable.end())
	      {
		it = ghostVertexTable.find(id);
		if (it==ghostVertexTable.end())
		  {
		    it = shadowVertexTable.find(id);
		    if (it==shadowVertexTable.end())
		      {
			ShadowVertex *sv;
			std::vector<Coord> pos(NDIM_W,0);
			LocalMesh::shadowVertexPool.pop(&sv);			
			st->getPosition(ITCell(0,id),&pos[0]);
			sv->set(&pos[0],
				LocalMesh::shadowVertexPool.getUsedCount()-1);

			shadowVertexTable.insert
			  (std::make_pair((unsigned long)id,
					  (void*)static_cast<Vertex*>(sv)));

			v[j]=static_cast<Vertex*>(sv);
			//printf("CREATING SHADOW VERTEX GEN = %ld\n",id);
		      } else v[j]=static_cast<Vertex*>(it->second);
		  } else v[j]=static_cast<Vertex*>(it->second);
	      } else v[j]=static_cast<Vertex*>(it->second);		    

	    v[j]->setSharedF();
	    v[j]->setGeneration(0,id);
	    v[j]->setGlobalIdentity(GlobalIdentity::EMPTY_RANK,id);
	  }
	shadow->setElements(v,projectedVolumeFunctor);
	shadow->setGlobalIdentity(GlobalIdentity::EMPTY_RANK,shadowId);
	shadow->setGeneration(0,shadowId);
	shadow->setLocalIndex(LocalMesh::shadowSimplexPool.getUsedCount()-1);
	shadowSimplexTable.insert
	  (std::make_pair(shadowId,(void*)static_cast<Simplex*>(shadow)));
      }
    else shadow = static_cast<ShadowSimplex*>(static_cast<Simplex*>(itS->second));

    // find the opposite vertex of the shadow simplex
    int difId=hlp::FindFirstDifference<ITCell,Simplex::NVERT>::find
      (&shadowVert[0],&neighborVert[0]);

    typename VT::iterator itV = vertexTable.find(shadowVert[difId].id());
    if (itV==vertexTable.end())
      itV = ghostVertexTable.find(shadowVert[difId].id());
    if (itV==ghostVertexTable.end()) 
      itV = shadowVertexTable.find(shadowVert[difId].id());
    if (itV==shadowVertexTable.end())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("That should not happen, there is a bug ;(\n");
	exit(-1);
      }

    return std::make_pair(shadow,static_cast<Vertex*>(itV->second));
  }

  // Assign a correct global identity to the shadow/ghost simplices and vertices
  // and send list of send/receive for shadow/ghost simplices.
  // This is used as a helper for initFromRegularGrid(...)
  template <class EA, class ST, class VT, class ES>
  void initExchangeCells(EA &exchangeArr, ST &simplexTable, VT &vertexTable, 
			  MpiCellDataExchangeT<Simplex,ES> &exchangeInfo)
  {    
    int myRank=mpiCom->rank();
    int nParts=mpiCom->size();

    exchangeInfo.send.resize(nParts);
    exchangeInfo.receive.resize(nParts);

    if (exchangeArr.size() == 0) return;

    for (int i=0;i<nParts;i++)
      {	    
	unsigned long len=exchangeArr.size();
	mpiCom->Bcast(&len,i,1);
	unsigned long totalLen=len*(1+Simplex::NVERT);
	    
	std::vector<GlobalIdentityValue> identityArr(totalLen);
	// locally get global identity as set from grid
	if (myRank==i)
	  {	
	    unsigned long dec=0;
	    for (unsigned long j=0;j<len;j++,dec+=(1+Simplex::NVERT))
	      {
		identityArr[dec] = exchangeArr[j]->getGlobalIdentity().get();
		for (int k=1;k<=Simplex::NVERT;k++)
		  identityArr[dec+k] = 
		    exchangeArr[j]->getVertex(k-1)->getGlobalIdentity().get();
	      }
	  }
	
	mpiCom->Bcast(identityArr,i);
	std::vector<unsigned char> remoteIndex(len*Simplex::NVERT,0);

	// retrieve corresponding remote simplex
	std::vector<unsigned long> nFound(nParts,0);
	if (myRank!=i)
	  {
	    unsigned long dec=0;
	    for (unsigned long j=0;j<len;j++,dec+=(1+Simplex::NVERT))
	      {
		typename ST::iterator it = 
		  simplexTable.find(GlobalIdentity(identityArr[dec]).id());

		if (it!=simplexTable.end())
		  {
		    Simplex *s=((Simplex*)it->second);
		    // this global identity is correct as it comes from a local simplex
		    identityArr[dec] = s->getGlobalIdentity(myRank).get();
		    exchangeInfo.send[i].push_back(s);
		    nFound[myRank]++;
		    // remoteIndex stores the order of vertices in the remote simplex
		    //int remoteIndex[Simplex::NVERT];
		    for (int k=0;k<Simplex::NVERT;k++)
		      {
			// retrieve the vertices as they were set in the local simplex
			typename VT::iterator it = 
			  vertexTable.find(GlobalIdentity(identityArr[dec+k+1]).id());

			if (it==vertexTable.end())
			  {
			    PRINT_SRC_INFO(LOG_ERROR);
			    glb::console->print<LOG_ERROR>
			      ("local and remote simplex do not mach.\n");
			    exit(-1);
			  }
			Vertex *v=(Vertex *)it->second;
			unsigned char *remoteId = &remoteIndex[j*Simplex::NVERT];
			for (remoteId[k]=0;remoteId[k]<Simplex::NVERT;remoteId[k]++) 
			  if (v==s->getVertex(remoteId[k])) break;

			// identity array contains the dummy global index if the vertex is 
			// shared or the actual true global index if it is a local vertex
			identityArr[dec+remoteId[k]+1] = v->getGlobalIdentity().get();
		      }
		  }
		else identityArr[dec] = 0;
	      }
	    if (nFound[myRank]) exchangeInfo.sendRank.push_back(i);
	  }
	else identityArr.assign(totalLen,0);
	    
	mpiCom->Allreduce_inplace(nFound,MPI_SUM); 
	mpiCom->Reduce_inplace(identityArr,i,MPI_SUM); 
	mpiCom->Reduce_inplace(remoteIndex,i,MPI_SUM); 

	// and set global IDs	   
	if (myRank==i)
	  {		
	    unsigned long dec=0;
	    for (unsigned long j=0;j<len;j++,dec+=(1+Simplex::NVERT))
	      {
		// set shadow simplex global ID
		exchangeArr[j]->setGlobalIdentity(identityArr[dec]);
		GlobalIdentity id(identityArr[dec]);
		exchangeInfo.receive[id.rank()].push_back(exchangeArr[j]);

		// order of vertices in the remote simplex		   
		unsigned char *remoteId = &remoteIndex[j*Simplex::NVERT];
		unsigned char swp[Simplex::NVERT];

		// reorder vertices in the shadow simplex
		for (int k=0;k<Simplex::NVERT;k++) swp[k]=k;
		for (int k=0;k<Simplex::NVERT;k++)
		  {
		    if (k==swp[remoteIndex[k]]) continue;
		    exchangeArr[j]->swapVertices(k,swp[remoteId[k]]);
		    std::swap(swp[k],swp[remoteId[k]]);
		  }		 		    
	      }
		
	    exchangeInfo.receiveRank.clear();
	    for (long j=0;j<nParts;j++) 
	      if (nFound[j]>0) exchangeInfo.receiveRank.push_back(j);		
	  }
      }    
    exchangeInfo.updateNCum();
  }
  
  unsigned long getRemoteNLocalCells(int rank, int type) const
  {
    return globalNLocalCells[type][rank];
  }

  unsigned long getGlobalNCells(int type) const
  {
    return globalNCellsCum[type].back();
  }
  
  double updateLoadImbalanceFactor()
  { 
    //unsigned long min=globalNLocalCells[NDIM][0];
    double avg=globalNLocalCells[NDIM][0];
    unsigned long max=globalNLocalCells[NDIM][0];
    
    for (long i=1;i<mpiCom->size();i++)
      {
	if (globalNLocalCells[NDIM][i]>max) max=globalNLocalCells[NDIM][i];
	avg+=globalNLocalCells[NDIM][i];
	//if (globalNLocalCells[NDIM][i]<min) min=globalNLocalCells[NDIM][i];
      }
    loadImbalanceFactor = (double(max)/avg)*mpiCom->size();
  
    return loadImbalanceFactor;    
  }
  
  void updateCellsCount()
  {
    std::vector<unsigned long> nCells = LocalMesh::getNCells();
    std::vector<unsigned long> nCellsTotal = LocalMesh::getNCellsTotal();
    std::vector<unsigned long> tmp;    
    
    int ct=0;
    for (int i=0;i<NDIM+1;i++) 
      if (nCells[i]!=0) ct++;

    tmp.assign(ct*(mpiCom->size()+1)*2,0);
    
    ct=0;
    for (int i=0;i<NDIM+1;i++) 
      {
	if (nCells[i]!=0)
	  {
	    long index = ct*(mpiCom->size()+1) + mpiCom->rank()+1;
	    tmp[index] = nCells[i];		   
	    ct++;
	  }
      }    

    mpiCom->Allreduce_inplace(tmp,MPI_SUM);

    ct=0;
    for (int i=0;i<NDIM+1;i++) 
      {
	if (nCells[i]!=0)
	  {
	    globalNCellsCum[i].assign(&tmp[ct*(mpiCom->size()+1)],
				      &tmp[(ct+1)*(mpiCom->size()+1)]);
	    globalNLocalCells[i].assign(&tmp[ct*(mpiCom->size()+1) + 1],
					&tmp[(ct+1)*(mpiCom->size()+1)]);
	    for (int j=0;j<mpiCom->size();j++)
	      globalNCellsCum[i][j+1]+=globalNCellsCum[i][j];
	    ct++;
	  }
      }
    
    updateLoadImbalanceFactor();
  }
  /*
  std::vector<long> getNShadowSimplicesSendCum()
  {
    std::vector<long> result(sendNeighborNodesRank.size()+1);
    result[0]=0;
    for (int i=0;i<sendNeighborNodesRank.size();i++) 
      result[i+1]=result[i]+shadowSimplicesSend[sendNeighborNodesRank[i]].size();
    return result;
  }
  */
  /*
  std::vector<long> getCumShadowSimplicesReceive()
  {
    std::vector<long> result(receiveNeighborNodesRank.size()+1);
    result[0]=0;
    for (int i=0;i<receiveNeighborNodesRank.size();i++) 
      result[i+1]=result[i]+shadowSimplicesReceive[receiveNeighborNodesRank[i]].size();
    return result;
    }
  */

  // WARNING: nCellsCum must be up to date (i.e. call updateCellsCount())
  GlobalIndex globalIdentityToGlobalIndex(GlobalIdentity identity, int type)
  {
    return globalNCellsCum[type][identity.rank()] + identity.id();
  }
  
  /*
  // FIX ME: fix xadj when there are boundaries !
  // FIX ME: add weights to account for refinment  
  void generateParmetisGraph(ParmetisParams &p)
  {    
    if (p.owner) p.freeData();
    
    p.owner=true;
    p.vtxdist=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*globalNCellsCum[NDIM].size());
    p.xadj=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*(LocalMesh::getNCells(NDIM)+1));
    p.adjncy=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*LocalMesh::getNCells(NDIM)*Simplex::NNEI);

    for (unsigned long i=0;i<globalNCellsCum[NDIM].size();i++)
      p.vtxdist[i]=globalNCellsCum[NDIM][i];

  
    #pragma omp parallel for
    for (unsigned long i=0;i<glb::num_omp_threads;i++)
      {	
	const simplexPtr_iterator it_end=LocalMesh::simplexEnd();
	for (simplexPtr_iterator it=LocalMesh::simplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {	    
	    Simplex *nei[Simplex::NNEI];
	    const int count = it->getNeighbors(nei);
	    const unsigned long j=it->getLocalIndex();

	    p.xadj[j]=j*Simplex::NNEI;	    
	    PartitionerIndex *adjncy=&p.adjncy[p.xadj[j]];	    

	    for (int k=0;k<count;k++) 
	      adjncy[k] = getGlobalIndex(nei[k]);	   
	  }
      }   
      p.xadj[LocalMesh::getNCells(NDIM)] = LocalMesh::getNCells(NDIM)*Simplex::NNEI;
    }
  */
};

/** \}*/
#include "../internal/namespace.footer"
#endif
