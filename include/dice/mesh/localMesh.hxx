#ifndef __LOCAL_MESH_HXX__
#define __LOCAL_MESH_HXX__

#include <vector>
#include <list>
#include <algorithm>
#include <functional>

#ifdef USE_GNU_PSORT
#include <parallel/algorithm>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../dice_globals.hxx"

#include "../tools/types/types.hxx"
#include "../tools/memory/memoryPool.hxx"
#include "../tools/memory/iterableMemoryPool.hxx"
#include "../tools/helpers/unionIterator.hxx"
#include "../tools/helpers/helpers.hxx"
#include "../tools/helpers/helpers_macros.hxx"
#include "../tools/helpers/loopUnroller.hxx"
#include "../tools/IO/myIO.hxx"
#include "../tools/IO/paramsParser.hxx"

#include "../mesh/simplex.hxx"
#include "../mesh/meshParams.hxx"
#include "../mesh/sharedTree.hxx"
#include "../mesh/cellData/cellDataFunctors.hxx"
#include "../mesh/cellData/cellDataFunctors_default.hxx"
#include "../mesh/meshIteratorsMacros.hxx"

#include "../geometry/boundaryType.hxx"
#include "../geometry/geometricProperties.hxx"
#include "../geometry/simplexVolume.hxx"
#include "../geometry/simplexUniformSampler.hxx"

#include "../IO/ndNetUnstructuredMesh.hxx"

#include "internal/getBall.hxx"


/**
 * @file 
 * @brief  A class used to define a local (non MPI) adaptive unstructured simplicial 
 * mesh embedded in a higher dimension space.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class LocalMeshT
 * \brief  A class used to define a local (non MPI) unstructured simplicial mesh embedded 
 * in a higher dimension space. This class is NOT designed to be used on its own, but 
 * should rather be inherited. See MeshT ...
 * In particular, this class defines iterators to the simplices/vertices and various 
 * helper functions used to access/modify the mesh.
 * \tparam T mesh traits 
 */

template <class T>
class LocalMeshT : public SharedTreeT< SimplexT<T> , T >
{  
public:
  typedef LocalMeshT<T> MyType;
  typedef SharedTreeT< SimplexT<T> , T > Tree;
  
  static const int NDIM                = T::NDIM;
  static const int NDIM_W              = T::NDIM_W;  
  static const int BOUNDARY_TYPE       = T::BOUNDARY_TYPE;
  static const int WORLD_BOUNDARY_TYPE = T::WORLD_BOUNDARY_TYPE;   

  static std::string classHeader() {return "local_mesh";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  typedef typename Tree::Element          Simplex;

  typedef typename Simplex::Vertex        Vertex; 
  typedef typename Simplex::Segment       Segment; 
  typedef typename Simplex::SegmentHandle SegmentHandle;
  typedef typename Simplex::Facet         Facet;
  typedef typename Simplex::FacetHandle   FacetHandle;

  typedef typename Simplex::Shadow        ShadowSimplex; 
  typedef typename Vertex::Shadow         ShadowVertex;
  typedef typename Simplex::Ghost         GhostSimplex;
  typedef typename Vertex::Ghost          GhostVertex;

  typedef IterableMemoryPoolT<Simplex,iteratorThreadModel::SortedBlocks> 
  SimplexPool;  
  typedef IterableMemoryPoolT<Vertex,iteratorThreadModel::SortedBlocks>        
  VertexPool;
  typedef IterableMemoryPoolT<GhostSimplex,iteratorThreadModel::SortedBlocks>  
  GhostSimplexPool;
  typedef IterableMemoryPoolT<GhostVertex,iteratorThreadModel::SortedBlocks>   
  GhostVertexPool;
  typedef IterableMemoryPoolT<ShadowSimplex,iteratorThreadModel::SortedBlocks> 
  ShadowSimplexPool; 
  typedef IterableMemoryPoolT<ShadowVertex,iteratorThreadModel::SortedBlocks>  
  ShadowVertexPool;

  typedef typename VertexPool::UnserializedPointerUpdater        UVPUpdater;
  typedef typename GhostVertexPool::UnserializedPointerUpdater   UGVPUpdater;
  typedef typename ShadowVertexPool::UnserializedPointerUpdater  USVPUpdater;
  typedef typename SimplexPool::UnserializedPointerUpdater       USPUpdater;
  typedef typename GhostSimplexPool::UnserializedPointerUpdater  UGSPUpdater;
  typedef typename ShadowSimplexPool::UnserializedPointerUpdater USSPUpdater;

  typedef typename Simplex::segment_circulator       segment_circulator;
  typedef typename Simplex::const_segment_circulator const_segment_circulator;
    
  // typedef typename Tree::Leaf Leaf;
  // typedef typename Tree::Root Root;
  typedef typename Tree::Node Node; 

  typedef typename T::GlobalIndex          GlobalIndex;
  typedef typename T::LocalIndex           LocalIndex;
  typedef typename Simplex::GlobalIdentity GlobalIdentity;
  typedef typename GlobalIdentity::Value   GlobalIdentityValue;

  typedef typename T::Coord Coord; 
  typedef T Traits;     

  typedef GeometricPropertiesT<Coord,NDIM,NDIM_W,T::BOUNDARY_TYPE,T::WORLD_BOUNDARY_TYPE> 
  GeometricProperties;
  
  typedef CellDataFunctorsT<MyType>                 CellDataFunctors;
  typedef typename CellDataFunctors::VertexFunctor  VertexFunctor;
  typedef typename CellDataFunctors::SimplexFunctor SimplexFunctor; 

  typedef std::pair<float,int> CheckRefineReturnType;

  // typedef SimplexVolumeT<NDIM,NDIM_W,Coord> VolumeFromCoords;
  // typedef SimplexVolumeT<NDIM,NDIM,Coord>   ProjectedVolumeFromCoords;
  typedef MeshParamsT<NDIM,NDIM_W,Coord>    Params;  
  /*
  static long getSimplexRefineBufferSize()
  {
    return T::template getSimplexRefineBufferSize<MyType,SegmentHandle,Vertex,Simplex>();
  }
  */
protected:
  typedef typename VertexPool::SortedPointerUpdater        SVPUpdater;
  typedef typename GhostVertexPool::SortedPointerUpdater   SGVPUpdater;
  typedef typename ShadowVertexPool::SortedPointerUpdater  SSVPUpdater;
  typedef typename SimplexPool::SortedPointerUpdater       SSPUpdater;
  typedef typename GhostSimplexPool::SortedPointerUpdater  SGSPUpdater;
  typedef typename ShadowSimplexPool::SortedPointerUpdater SSSPUpdater;

  Params params;

  LocalMeshT(MpiCommunication *com):    
    Tree(com),     
    mpiCom(com),
    simplexPool("Simplex"),    
    vertexPool("Vertex"),
    ghostSimplexPool("GhostSimplex"),
    ghostVertexPool("GhostVertex"),
    shadowSimplexPool("ShadowSimplex"),
    shadowVertexPool("ShadowVertex"),
    cellDataFunctors(this),
    geometry(NULL),
    initialized(false)
  {
    
  }
  
  /*
  template <class R>
  LocalMeshT(const MeshParams &defParams, MpiCommunication *com, 
	     ParamsParser *parser, R *reader):    
    Tree(com,parser,reader),     
    mpiCom(com),
    simplexPool("Simplex"),    
    vertexPool("Vertex"),
    ghostSimplexPool("GhostSimplex"),
    ghostVertexPool("GhostVertex"),
    shadowSimplexPool("ShadowSimplex"),
    shadowVertexPool("ShadowVertex"),
    cellDataFunctors(this),
    initialized(false)
  {
  
  }
*/

  /*

  template <class R>
  LocalMeshT(MpiCommunication *com, 
	     ParamsParser *parser, R *reader):    
    Tree(com,parser,reader),     
    mpiCom(com),
    simplexPool("Simplex"),    
    vertexPool("Vertex"),
    ghostSimplexPool("GhostSimplex"),
    ghostVertexPool("GhostVertex"),
    shadowSimplexPool("ShadowSimplex"),
    shadowVertexPool("ShadowVertex"),
    cellDataFunctors(this),
    initialized(false)
  {
  
  }
  */
  virtual ~LocalMeshT()
  {
    delete geometry;  
  }

  void defrag()
  {
    UVPUpdater vpu;
    UGVPUpdater gvpu;
    USVPUpdater svpu;

    USPUpdater spu;
    UGSPUpdater gspu;
    USSPUpdater sspu;

    bool swap=false;
    spu  = simplexPool.defrag();
    vpu  = vertexPool.defrag();
    gspu = ghostSimplexPool.defrag();
    gvpu = ghostVertexPool.defrag();
    sspu = shadowSimplexPool.defrag();
    svpu = shadowVertexPool.defrag();

#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const simplexPtr_iterator it_end=simplexEnd();
	for (simplexPtr_iterator it=simplexBegin(i,glb::num_omp_threads);
	     it!=it_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);
	
	const ghostSimplexPtr_iterator itg_end=ghostSimplexEnd();
	for (ghostSimplexPtr_iterator it=ghostSimplexBegin(i,glb::num_omp_threads);
	     it!=itg_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);
	
	const shadowSimplexPtr_iterator its_end=shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=shadowSimplexBegin(i,glb::num_omp_threads);
	     it!=its_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);	 
      } 

    Tree::defrag(spu);
  }

  void construct()
  {    
    int flags = cellDataFunctors::F_NO_FLAG;
    // register mock data for simplices and vertices   
    if (NDIM_W>NDIM) cellDataFunctors.template insert<cellDataFunctors::ExtraDimsT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SignedVolumeT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SignedProjectedVolumeT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::VertexGenerationT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::DomainIndexT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::LevelT>(flags);

    // We do not dump these functors to files unless we are debugging
    if (!glb::debug) flags |= cellDataFunctors::F_SKIP_ON_FILE_DUMP;
    
    cellDataFunctors.template insert<cellDataFunctors::VertexFlagsT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SimplexFlagsT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::VertexLocalIndexT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SimplexLocalIndexT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::VertexGlobalIndexT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SimplexGlobalIndexT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::VertexRankT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::SimplexRankT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::WeightT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::DepthT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::ProjectedAnisotropyT>(flags);
    cellDataFunctors.template insert<cellDataFunctors::ProjectedMinHeightT>(flags);

    // volumeFunctor = cellDataFunctors.getSimplexData("volume");
    // projectedVolumeFunctor = cellDataFunctors.getSimplexData("projectedVolume");
  }  
  
  void build(const Params &meshParams)
  {
    params = meshParams;
           
    simplexPool.setAllocFactor(params.allocFactor);
    vertexPool.setAllocFactor(params.allocFactor);
    ghostSimplexPool.setAllocFactor(params.allocFactor);
    ghostVertexPool.setAllocFactor(params.allocFactor);
    shadowSimplexPool.setAllocFactor(params.allocFactor);
    shadowVertexPool.setAllocFactor(params.allocFactor);
    Tree::setAllocFactor(params.allocFactor);

    if (geometry != NULL) delete geometry;
    geometry = new GeometricProperties(&params.x0[0],&params.delta[0]);	
  }

  template <class R>
  void build(const Params &meshParams, R *reader)
  {
    UVPUpdater vpu;
    UGVPUpdater gvpu;
    USVPUpdater svpu;

    USPUpdater spu;
    UGSPUpdater gspu;
    USSPUpdater sspu;

    build(meshParams,reader,vpu,gvpu,svpu,spu,gspu,sspu);    
  }

  template <class R>
  void build(const Params &meshParams, R *reader,
	     UVPUpdater &vpu,
	     UGVPUpdater &gvpu,
	     USVPUpdater &svpu,
	     USPUpdater &spu,
	     UGSPUpdater &gspu,
	     USSPUpdater &sspu)
  {  
    build(meshParams);
    if (reader == NULL) return;    
    
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);
    
    // restore the memory pools here   
    glb::console->printFlush<LOG_INFO>("Unserializing vertices pools ... ");
    vpu = vertexPool.unSerialize(reader);   
    gvpu = ghostVertexPool.unSerialize(reader);
    svpu = shadowVertexPool.unSerialize(reader);
    glb::console->printFlush<LOG_INFO>("done.\n");

    glb::console->printFlush<LOG_INFO>("Unserializing simplices pools ... ");
    spu = simplexPool.unSerialize(reader);    
    gspu = ghostSimplexPool.unSerialize(reader);
    sspu = shadowSimplexPool.unSerialize(reader);
    glb::console->printFlush<LOG_INFO>("done.\n");
    
    glb::console->printFlush<LOG_INFO>("Updating cell pointers ... ");
    bool swap = reader->getNeedSwap();
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const simplexPtr_iterator it_end=simplexEnd();
	for (simplexPtr_iterator it=simplexBegin(i,glb::num_omp_threads);
	     it!=it_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);
	
	const ghostSimplexPtr_iterator itg_end=ghostSimplexEnd();
	for (ghostSimplexPtr_iterator it=ghostSimplexBegin(i,glb::num_omp_threads);
	     it!=itg_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);
	 
	const shadowSimplexPtr_iterator its_end=shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=shadowSimplexBegin(i,glb::num_omp_threads);
	     it!=its_end;++it)
	  (*it)->updateAfterUnserialized(*vpu,*gvpu,*svpu,*spu,*gspu,*sspu,swap);	 
      } 
    glb::console->printFlush<LOG_INFO>("done.\n");

    Tree::build(reader,spu);
  }

public:

  const SimplexFunctor *getSimplexFunctorPtr(const std::string &name)
  {
    SimplexFunctor *d=cellDataFunctors.getSimplexFunctor(name);
    return d;
  }

  const VertexFunctor *getVertexFunctorPtr(const std::string &name)
  {
    VertexFunctor *d=cellDataFunctors.getVertexFunctor(name);
    return d;
  }

  const SimplexFunctor *getSimplexFunctorPtr(int i)
  {
    SimplexFunctor *d=cellDataFunctors.getSimplexFunctor(i);
    return d;
  }

  const VertexFunctor *getVertexFunctorPtr(int i)
  {
    VertexFunctor *d=cellDataFunctors.getVertexFunctor(i);
    return d;
  }

  int getNSimplexFunctor() const
  {
    return cellDataFunctors.getSimplexFunctorCount();
  }

  int getNVertexFunctor() const
  {
    return cellDataFunctors.getVertexFunctorCount();
  }

  template <template <class CFM> class CFT>
  bool addCellDataFunctor(int flags=cellDataFunctors::F_NO_FLAG, 
			  bool replaceIfExists=false)
  {
    return cellDataFunctors.template insert<CFT>(flags,replaceIfExists);
  }

  GeometricProperties *getGeometry() const
  {
    return geometry;
  }

  std::vector<unsigned long> getNCells()
  {
    std::vector<unsigned long> result(NDIM+1,0);
    result[0]=vertexPool.getUsedCount();
    //result[1]=segmentPool.getUsedCount();
    result[NDIM]=simplexPool.getUsedCount();    
    return result;
  }

  std::vector<unsigned long> getNCellsTotal()
  {
    std::vector<unsigned long> result(NDIM+1,0);
    result[0]=vertexPool.getUsedCount()+
      shadowVertexPool.getUsedCount()+
      ghostVertexPool.getUsedCount();
    
    result[NDIM]=simplexPool.getUsedCount()+
      shadowSimplexPool.getUsedCount()+
      ghostSimplexPool.getUsedCount();

    return result;
  }

  unsigned long getNCells(int i)
  {
    if (i==0) 
      return vertexPool.getUsedCount();
    
    if (i==NDIM) 
      return simplexPool.getUsedCount();    
    return 0;
  }

  unsigned long getNCellsTotal(int i)
  {
    if (i==0) 
      return vertexPool.getUsedCount()+
	shadowVertexPool.getUsedCount()+
	ghostVertexPool.getUsedCount();

    if (i==NDIM) 
      return simplexPool.getUsedCount()+
	shadowSimplexPool.getUsedCount()+
	ghostSimplexPool.getUsedCount();   

    return 0;
  }
  /*
  std::vector<unsigned long> getGlobalNCells()
  {
    return globalNCells;
  }

  unsigned long getGlobalNCells(int i)
  {
    if (i<=NDIM) return globalNCells[i];
    return 0;
  }
  */
  unsigned long getNSimplices()
  {
    return simplexPool.getUsedCount();
  }

  unsigned long getNVertices()
  {
    return vertexPool.getUsedCount();
  }

  unsigned long getNShadowSimplices()
  {
    return shadowSimplexPool.getUsedCount();
  }

  unsigned long getNShadowVertices()
  {
    return shadowVertexPool.getUsedCount();
  }

  unsigned long getNGhostSimplices()
  {
    return ghostSimplexPool.getUsedCount();
  }

  unsigned long getNGhostVertices()
  {
    return ghostVertexPool.getUsedCount();
  }
  
  std::vector<Vertex*> getVerticesArray()
  {
    std::vector<Vertex*> result(vertexPool.getUsedCount());
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const vertexPtr_iterator it_end=vertexEnd();
	for (vertexPtr_iterator it=vertexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }
  
  std::vector<Vertex*> getSharedVerticesArray()
  {
    std::vector<Vertex*> result;    
    long ns =shadowSimplexPool.getUsedCount()+ghostSimplexPool.getUsedCount();
    if (ns == 0) return result;
    result.reserve(ns);

    // only few regular vertices are shared so this is fast
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const vertexPtr_iterator it_end=vertexEnd();
	for (vertexPtr_iterator it=vertexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    if (it->isShared())
	      {
#pragma omp critical
		result.push_back(*it);
	      }
	  }
      }
    
    // but all ghost/shadow vertices are shared so no use for openMP
    const ghostVertexPtr_iterator git_end=ghostVertexEnd();
    for (ghostVertexPtr_iterator git=ghostVertexBegin();git!=git_end;++git)
      result.push_back((Vertex*)(*git));
   
    const shadowVertexPtr_iterator sit_end=shadowVertexEnd();
    for (shadowVertexPtr_iterator sit=shadowVertexBegin();sit!=sit_end;++sit)
      result.push_back((Vertex*)(*sit));
     
    return result;    
  }
  
  std::vector<Simplex*> getSimplicesArray()
  {
    std::vector<Simplex*> result(simplexPool.getUsedCount());
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const simplexPtr_iterator it_end=simplexEnd();
	for (simplexPtr_iterator it=simplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }
  
  template <class Functor>
  void visitSimplices(const Functor &f, 
		      bool visitLocals, bool visitGhosts, bool visitShadows,
		      int nThreads=glb::num_omp_threads)
  {
    if (nThreads<1) nThreads=omp_get_max_threads();
    if (visitLocals)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {
	    int th=omp_get_thread_num();
	    const simplexPtr_iterator it_end=simplexEnd();
	    for (simplexPtr_iterator it=simplexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
    if (visitGhosts)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {
	    int th=omp_get_thread_num();
	    const ghostSimplexPtr_iterator it_end=ghostSimplexEnd();
	    for (ghostSimplexPtr_iterator it=ghostSimplexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
    if (visitShadows)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {
	    int th=omp_get_thread_num();
	    const shadowSimplexPtr_iterator it_end=shadowSimplexEnd();
	    for (shadowSimplexPtr_iterator it=shadowSimplexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
  }
  
  template <class Functor>
  void visitVertices(Functor &f, 
		     bool visitLocals, bool visitGhosts, bool visitShadows,
		     int nThreads=glb::num_omp_threads)
  {
    if (nThreads<1) nThreads=omp_get_max_threads();
    if (visitLocals)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {	   
	    int th=omp_get_thread_num();
	    const vertexPtr_iterator it_end=vertexEnd();
	    for (vertexPtr_iterator it=vertexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
    if (visitGhosts)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {
	    int th=omp_get_thread_num();
	    const ghostVertexPtr_iterator it_end=ghostVertexEnd();
	    for (ghostVertexPtr_iterator it=ghostVertexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
    if (visitShadows)
      {
#pragma omp parallel for num_threads(nThreads)
	for (long i=0;i<nThreads;i++)
	  {
	    int th=omp_get_thread_num();
	    const shadowVertexPtr_iterator it_end=shadowVertexEnd();
	    for (shadowVertexPtr_iterator it=shadowVertexBegin(i,nThreads);
		 it!=it_end;++it)
	      {
		f(*it,th);
	      }
	  }
      }
  }
  
  /*
  virtual void dumpGhostsToAscii(const std::string &fname_, bool completeFName=false)
  {
    if (getNGhostSimplices()==0) return;

    std::string fname(fname_);
    if (completeFName)
      {
	char tmp[255];
	sprintf(tmp,"ghost_%s_%4.4d.dat",fname.c_str(),mpiCom->rank());
	fname = std::string(tmp);
      }

    FILE *f=fopen(fname.c_str(),"w");
    const ghostSimplexPtr_iterator it_end=ghostSimplexEnd();
    for (ghostSimplexPtr_iterator it=ghostSimplexBegin();it!=it_end;++it)
      {
	Simplex *s=(*it);
	Coord *pos[Simplex::NVERT];
	s->getVerticesCoordsPtr(pos);
	for (int i=0;i<Simplex::NVERT;i++)
	  {
	    for (int j=0;j<NDIM_W;j++)
	      fprintf(f,"%lg ",pos[i][j]);	    
	  }
	fprintf(f,"%ld \n",(long)it->getLocalIndex());
      }
    fclose(f);
  }

  virtual void dumpShadowsToAscii(const std::string &fname_, bool completeFName=false)
  {
    if (getNShadowSimplices()==0) return;

    std::string fname(fname_);
    if (completeFName)
      {
	char tmp[255];
	sprintf(tmp,"shadow_%s_%4.4d.dat",fname.c_str(),mpiCom->rank());
	fname = std::string(tmp);
      }  

    FILE *f=fopen(fname.c_str(),"w");
    const shadowSimplexPtr_iterator it_end=shadowSimplexEnd();
    for (shadowSimplexPtr_iterator it=shadowSimplexBegin();it!=it_end;++it)
      {
	Simplex *s=(*it);
	Coord *pos[Simplex::NVERT];
	s->getVerticesCoordsPtr(pos);
	for (int i=0;i<Simplex::NVERT;i++)
	  {
	    for (int j=0;j<NDIM_W;j++)
	      fprintf(f,"%lg ",pos[i][j]);	    
	  }
	fprintf(f,"%ld \n",(long)it->getLocalIndex());
      }
    fclose(f);
  }
  */

  virtual void dumpToNDnetwork_DBG(const std::string &fname_, bool completeFName=false, 
				   bool dumpAll=false, bool withNeighbors=false)
  {
    std::vector<unsigned long> nCells = (dumpAll)?getNCellsTotal():getNCells();    
    char comment[80];
    std::string fname(fname_);
   
    if ((completeFName)&&(mpiCom->size()>1))
      {
	char tmp[255];
	sprintf(tmp,"%s_%4.4d.NDnet",fname.c_str(),mpiCom->rank());
	fname = std::string(tmp);
      }  
    else
      {
	char tmp[255];
	sprintf(tmp,"%s.NDnet",fname.c_str());
	fname = std::string(tmp);
      }

    sprintf(comment,"Local mesh %d / %d",(int)mpiCom->rank(),(int)mpiCom->size());
    
    IO::NDnetwork net(NDIM, NDIM_W,&params.x0[0],&params.delta[0],&nCells[0],
		      comment,BOUNDARY_TYPE==BoundaryType::PERIODIC);
    std::vector<Vertex*> vArr=getVerticesArray();
    std::vector<GhostVertex*> gvArr;
    std::vector<ShadowVertex*> svArr;
    if (dumpAll)
      {
	gvArr=getGhostVerticesArray();
	svArr=getShadowVerticesArray();
      }

    if (mpiCom->size()==1)
      glb::console->printFlush<LOG_STD>
	("Dumping local %dD mesh to %dDnetwork file '%s' ... ",NDIM,NDIM_W,fname.c_str());
    else
      glb::console->printFlush<LOG_STD>
	("Dumping local %dD mesh to %dDnetwork files '%s' (%d files) ... ",
	 NDIM,NDIM_W,fname.c_str(),mpiCom->size());

    FILE *f=fopen(fname.c_str(),"w");
    if (f==NULL)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("Could not open file '%s' for writing.\n",fname.c_str());
	exit(-1);
      }

     typename myIO::BinaryWriterT<> bWrite(f);
    // header
    net.writeHeader(f);
    
    // vertex coordinates
    unsigned int jj=net.ndims*nCells[0];
    fwrite(&jj,sizeof(unsigned int),1,f);
    
    for (unsigned long i=0;i<vArr.size();i++)
      {
	bWrite.template writeAs<NDNET_FLOAT>(vArr[i]->coords,NDIM_W);	
      }
    if (dumpAll)
      {
	for (unsigned long i=0;i<gvArr.size();++i)
	  bWrite.template writeAs<NDNET_FLOAT>(gvArr[i]->coords,NDIM_W);
	for (unsigned long i=0;i<svArr.size();++i)
	  bWrite.template writeAs<NDNET_FLOAT>(svArr[i]->coords,NDIM_W);
      }
    bWrite.flush();
    fwrite(&jj,sizeof(unsigned int),1,f);
  
    // cell count
    jj=(1+net.ndims)*sizeof(NDNET_UINT);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(net.nfaces,sizeof(NDNET_UINT),((size_t)net.ndims+1),f);
    fwrite(&jj,sizeof(unsigned int),1,f);

    // defined cells
    jj=(1+net.ndims)*sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(net.haveVertexFromFace,sizeof(int),((size_t)net.ndims+1),f);
    fwrite(&jj,sizeof(unsigned int),1,f);

    // simplices
    std::vector<Simplex*> sArr;  
    std::vector<GhostSimplex*> gsArr;  
    std::vector<ShadowSimplex*> ssArr;  
    const simplexPtr_iterator sim_end=simplexEnd();
    jj=sizeof(NDNET_UINT)*((size_t)(NDIM+1)*nCells[NDIM]);
    fwrite(&jj,sizeof(unsigned int),1,f);

    if (withNeighbors)
      {
	sArr=getSimplicesArray();	
	for (unsigned long j=0;j<sArr.size();++j)
	  {	
	     //NDNET_UINT tmp[(NDIM+1)];
	    typename Vertex::LocalIndex tmp[Simplex::NVERT];
	    sArr[j]->getVerticesLocalIndex(tmp);
	    bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
	  }       
      }
    else
      {	
	for (simplexPtr_iterator it=simplexBegin();it!=sim_end;++it)
	  {
	     //NDNET_UINT tmp[(NDIM+1)];
	    typename Vertex::LocalIndex tmp[Simplex::NVERT];
	    it->getVerticesLocalIndex(tmp);	
	    bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
	  }
      }

    if (dumpAll)
      {
	gsArr=getGhostSimplicesArray();
	ssArr=getShadowSimplicesArray();
	long delta=getNVertices();
	for (unsigned long j=0;j<gsArr.size();++j)
	  {
	    //NDNET_UINT tmp[(NDIM+1)];
	    typename Vertex::LocalIndex tmp[Simplex::NVERT];
	    gsArr[j]->getVerticesLocalIndex(tmp);
	    for (int k=0;k<NDIM+1;k++) 
	      {
		if (gsArr[j]->getVertex(k)->isGhost()) tmp[k]+=delta;
	      }
	    bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
	  }   	
	for (unsigned long j=0;j<ssArr.size();++j)
	  {
	    //NDNET_UINT tmp[(NDIM+1)];
	    typename Vertex::LocalIndex tmp[Simplex::NVERT];
	    ssArr[j]->getVerticesLocalIndex(tmp);
	    for (int k=0;k<NDIM+1;k++) 
	      {
		if (ssArr[j]->getVertex(k)->isGhost()) tmp[k]+=delta;
		else if (ssArr[j]->getVertex(k)->isShadow()) tmp[k]+=delta+gvArr.size();
	      }
	    bWrite.template writeAs<NDNET_UINT>(tmp,(NDIM+1));	
	  }       
      }
    
    bWrite.flush();
    fwrite(&jj,sizeof(unsigned int),1,f);    

    // junk
    jj=(1+net.ndims)*sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(net.haveFaceFromVertex,sizeof(int),((size_t)net.ndims+1),f);
    fwrite(&jj,sizeof(unsigned int),1,f);
        
    if (withNeighbors) net.haveFaceFromFace[NDIM][NDIM]=1;  

    jj=(1+net.ndims)*(1+net.ndims)*sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    for (long i=0;i<net.ndims+1;i++)
      fwrite(net.haveFaceFromFace[i],sizeof(int),((size_t)net.ndims+1),f);
    fwrite(&jj,sizeof(unsigned int),1,f);

    if (withNeighbors) 
      {	
	int i=NDIM;	
	std::vector<NDNET_IDCUMT> nn(nCells[i]+1);
	nn[0]=0;
	for (unsigned long j=1;j<nn.size();++j)
	  nn[j]=nn[j-1]+Simplex::NNEI;
	jj=sizeof(NDNET_IDCUMT)*nn.size();
	fwrite(&jj,sizeof(unsigned int),1,f);
	fwrite(&nn[0],sizeof(NDNET_IDCUMT),((size_t)net.nfaces[i]+1),f);
	fwrite(&jj,sizeof(unsigned int),1,f);
	nn.clear();
	
	jj=sizeof(NDNET_UINT)*(nCells[i]*Simplex::NNEI);
	fwrite(&jj,sizeof(unsigned int),1,f);
	const long myRank=mpiCom->rank();
	for (unsigned long j=0;j<sArr.size();++j)
	  {
	    NDNET_UINT tmp[Simplex::NNEI];
	    Simplex *cur=sArr[j];
	    for (int k=0;k<Simplex::NNEI;++k)
	      {
		//Simplex *cur=sArr[j];
		Simplex *nei=cur->getNeighbor(k);
		if (nei==NULL)
		  tmp[k]=j;
		else if (myRank!=nei->getGlobalIdentity(myRank).rank())
		  tmp[k]=j;
		else
		  tmp[k]=nei->getLocalIndex();
	      }
	    bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
	  }
	if (dumpAll)
	  {
	    long delta[3];
	    delta[0]=0;
	    delta[1]=sArr.size();
	    delta[2]=delta[1]+gsArr.size();
	    for (unsigned long j=0;j<gsArr.size();++j)
	      {
		NDNET_UINT tmp[Simplex::NNEI];
		for (int k=0;k<Simplex::NNEI;++k)
		  {
		    Simplex *nei=gsArr[j]->getNeighbor(k);
		    if (nei==NULL)
		      tmp[k]=delta[1]+j;
		    else 
		      {
			tmp[k]=nei->getLocalIndex();
			if (nei->isGhost())
			  tmp[k]+=delta[1];
			else if (nei->isShadow())
			  tmp[k]+=delta[2]; 	    
		      }
		  }
		bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
	      }

	    for (unsigned long j=0;j<ssArr.size();++j)
	      {
		NDNET_UINT tmp[Simplex::NNEI];
		for (int k=0;k<Simplex::NNEI;++k)
		  {
		    Simplex *nei=ssArr[j]->getNeighbor(k);
		    if (nei==NULL)
		      tmp[k]=delta[2]+j;
		    else 
		      {
			tmp[k]=nei->getLocalIndex();
		    
			if (nei->isGhost())
			  tmp[k]+=delta[1];
			else if (nei->isShadow())
			  tmp[k]+=delta[2]; 	
		      }
		  }
		bWrite.template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
	      }
	  }
	bWrite.flush();
	fwrite(&jj,sizeof(unsigned int),1,f);
      }

    bWrite.flush();
    net.haveVFlags=1;
    jj=sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(&net.haveVFlags,sizeof(int),1,f);
    fwrite(&jj,sizeof(unsigned int),1,f);
      
    // vertex flags
    jj=sizeof(unsigned char)*nCells[0];
    fwrite(&jj,sizeof(unsigned int),1,f);
    for (unsigned long i=0;i<vArr.size();i++)
      bWrite.template writeAs<unsigned char>(&(vArr[i]->flags));
    if (dumpAll)
      {
	for (unsigned long i=0;i<gvArr.size();i++)
	  bWrite.template writeAs<unsigned char>(&(gvArr[i]->flags));
	for (unsigned long i=0;i<svArr.size();i++)
	  bWrite.template writeAs<unsigned char>(&(svArr[i]->flags));
      }
    bWrite.flush();
    fwrite(&jj,sizeof(unsigned int),1,f);
    
    // cell flags

    net.haveFFlags[NDIM]=1;
    jj=sizeof(int)*(net.ndims+1);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(net.haveFFlags,sizeof(int),(net.ndims+1),f);
    fwrite(&jj,sizeof(unsigned int),1,f);
  
    //simplex flags
    jj=sizeof(unsigned char)*nCells[NDIM];
    fwrite(&jj,sizeof(unsigned int),1,f);
    if (withNeighbors)
      {
	for (unsigned long i=0;i<sArr.size();++i)
	  bWrite.template writeAs<unsigned char>(&(sArr[i]->flags));	
      }
    else
      {
	for (simplexPtr_iterator it=simplexBegin();it!=sim_end;++it)
	  bWrite.template writeAs<unsigned char>(&(it->flags));	
      }
    if (dumpAll)
      {
	for (unsigned long i=0;i<gsArr.size();i++)
	  bWrite.template writeAs<unsigned char>(&(gsArr[i]->flags));
	for (unsigned long i=0;i<ssArr.size();i++)
	  bWrite.template writeAs<unsigned char>(&(ssArr[i]->flags));
      }
    bWrite.flush();
    fwrite(&jj,sizeof(unsigned int),1,f);
    
    // data
    int nVertexData=cellDataFunctors.getVertexFunctorCount();
    int nSimplexData=cellDataFunctors.getSimplexFunctorCount();
    
    int nData=0;
    for (int i=0;i<nVertexData;i++) 
      {
	const auto *vf=getVertexFunctorPtr(i);
	if (!(vf->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP))
	  nData+=vf->getSize();      
      }
    for (int i=0;i<nSimplexData;i++) 
      {
	const auto *sf=getSimplexFunctorPtr(i);
	if (!(sf->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP))
	  nData+=sf->getSize();
      }
    /*
    for (int i=0;i<nVertexData;i++) 
      nData+=cellDataFunctors.getVertexFunctor(i)->getSize();
    for (int i=0;i<nSimplexData;i++) 
      nData+=cellDataFunctors.getSimplexFunctor(i)->getSize();
    */
    jj=sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(&nData,sizeof(int),1,f);
    fwrite(&jj,sizeof(unsigned int),1,f);

    if (nData)
      {	
	int type=0;
	int count=0;
	char name[255];
	int nData = cellDataFunctors.getVertexFunctorCount();

	for (int i=0;i<nData;i++)
	  {
	    VertexFunctor *gvd=cellDataFunctors.getVertexFunctor(i);
	    count=gvd->getSize();

	    if (gvd->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP)
		continue;

	    for (int ct=0;ct<count;ct++)
	      {//printf(name,"%s_%2.2d\n",gsd->getName().c_str(),ct);
		if (count>1)
		  sprintf(name,"%s_%2.2d",gvd->getName().c_str(),ct);
		else
		  strcpy(name,gvd->getName().c_str());
		
		jj=sizeof(int)+255*sizeof(char);
		fwrite(&jj,sizeof(unsigned int),1,f);
		fwrite(&type,sizeof(int),1,f);
		fwrite(name,sizeof(char)*255,1,f);
		fwrite(&jj,sizeof(unsigned int),1,f);

		jj=sizeof(double)*nCells[0];
		fwrite(&jj,sizeof(unsigned int),1,f);
		
		for (unsigned long j=0;j<vArr.size();j++)
		  {
		    double tmp=gvd->get(vArr[j],ct);
		    bWrite.template writeAs<double>(&tmp);
		  }
	
		if (dumpAll)
		  {
		    for (unsigned long j=0;j<gvArr.size();j++)
		      {
			double tmp=gvd->get(gvArr[j],ct);
			bWrite.template writeAs<double>(&tmp);
		      }
		    for (unsigned long j=0;j<svArr.size();j++)
		      {
			double tmp=gvd->get(svArr[j],ct);
			bWrite.template writeAs<double>(&tmp);
		      }		    
		  }
		bWrite.flush();
		fwrite(&jj,sizeof(unsigned int),1,f);
	      }		  
	  }
	
	bWrite.flush();
	nData = cellDataFunctors.getSimplexFunctorCount();
	type=NDIM;
	
	for (int i=0;i<nData;i++)
	  {
	    SimplexFunctor *gsd=cellDataFunctors.getSimplexFunctor(i);
	    count=gsd->getSize();

	    if (gsd->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP)
		continue;

	    for (int ct=0;ct<count;ct++)
	      {//printf(name,"%s_%2.2d\n",gsd->getName().c_str(),ct);
		if (count>1)
		  sprintf(name,"%s_%2.2d",gsd->getName().c_str(),ct);
		else 
		  strcpy(name,gsd->getName().c_str());
		
		jj=sizeof(int)+255*sizeof(char);
		fwrite(&jj,sizeof(unsigned int),1,f);
		fwrite(&type,sizeof(int),1,f);
		fwrite(name,sizeof(char)*255,1,f);
		fwrite(&jj,sizeof(unsigned int),1,f);

		//printf("WRITING field %s (%ld ele).\n",name,);

		jj=sizeof(double)*nCells[NDIM-1];
		fwrite(&jj,sizeof(unsigned int),1,f);
		if (withNeighbors)
		  {
		    for (unsigned long j=0;j<sArr.size();++j)
		      {
			double tmp=gsd->get(sArr[j],ct);
			bWrite.template writeAs<double>(&tmp);	
		      }		  
		  }
		else
		  {
		    for (simplexPtr_iterator it=simplexBegin();it!=sim_end;++it)
		      {
			double tmp=gsd->get(*it,ct);
			bWrite.template writeAs<double>(&tmp);		
		      }		    
		  }
		if (dumpAll)
		  {
		    for (unsigned long j=0;j<gsArr.size();j++)
		      {
			double tmp=gsd->get(gsArr[j],ct);
			bWrite.template writeAs<double>(&tmp);
		      }
		    for (unsigned long j=0;j<ssArr.size();j++)
		      {
			double tmp=gsd->get(ssArr[j],ct);
			bWrite.template writeAs<double>(&tmp);
		      }
		  }
		bWrite.flush();
		fwrite(&jj,sizeof(unsigned int),1,f);
	      }		  
	  }
      }
    bWrite.flush();
    //supdata
    jj=sizeof(int);
    fwrite(&jj,sizeof(unsigned int),1,f);
    fwrite(&net.nsupData,sizeof(int),1,f);
    fwrite(&jj,sizeof(unsigned int),1,f);

    fclose(f);
    glb::console->print<LOG_STD>("done.\n");
  }
  /*
  double distance(const Vertex *a,const Vertex *b) const
  {
    return geometry->distance(a->getCoordsConstPtr(),
			      b->getCoordsConstPtr());
  }

  double distance2(const Vertex *a,const Vertex *b) const
  {
    return geometry->distance2(a->getCoordsConstPtr(),
			       b->getCoordsConstPtr());
  }
  
  double length(SegmentHandle &s) const
  {
    return geometry->distance(s->getVertex(0)->getCoordsConstPtr(),
			      s->getVertex(1)->getCoordsConstPtr());
  }

  double length2(SegmentHandle &s) const
  {
    return geometry->distance2(s->getVertex(0)->getCoordsConstPtr(),
			       s->getVertex(1)->getCoordsConstPtr());
  }

  void checkBoundary(Vertex *v) const
  {
    geometry->checkBoundary(v->getCoordsPtr());
  }

  template <class OT>
  void midPointCoords(const Vertex *a,const Vertex *b, OT* out) const
  {
    geometry->midPointCoords(a->getCoordsConstPtr(),
			     b->getCoordsConstPtr(),
			     out); 
  }

  template <class OT>
  void midPointCoords(SegmentHandle &seg, OT* out) const
  {
    geometry->midPointCoords(seg->getVertex(0)->getCoordsConstPtr(),
			     seg->getVertex(1)->getCoordsConstPtr(),
			     out); 
  }
  */

  template <class S, class TT, int D=S::NDIM_W>
  void getBaseVectors(const S *s,TT (&vec)[S::NVERT-1][D]) const
  {
    const Coord *p[S::NVERT];
    s->getVerticesCoordsConstPtr(p);
    geometry->template getBaseVectors<Coord,TT,S::NVERT-1,D>(p,vec);    
  }

  template <class S,class TT>
  void getProjectedBaseVectors(const S *s, TT (&vec)[S::NVERT-1][NDIM]) const
  {
    getBaseVectors<S,TT,S::NDIM>(s,vec);    
  }

  void computeLocalBoundingBox(double xMin[NDIM], double xMax[NDIM], 
			       int nThreads=glb::num_omp_threads)
  {
    if (nThreads<1) nThreads=omp_get_max_threads();

    for (int i=0;i<NDIM;++i) 
      {
	xMin[i]=std::numeric_limits<double>::max();
	xMax[i]=-std::numeric_limits<double>::max();
      }

    const Coord *refCoord = vertexBegin()->getCoordsConstPtr();
   
#pragma omp parallel for num_threads(nThreads)
    for (long th=0;th<nThreads;th++)
      {
	double min[NDIM];
	double max[NDIM];
	Coord ref[NDIM];

	for (int i=0;i<NDIM;++i) 
	  {
	    ref[i]=refCoord[i];
	    min[i]=std::numeric_limits<double>::max();
	    max[i]=-std::numeric_limits<double>::max();
	  }
	
	const vertexPtr_LGS_iterator itv_end=vertexLGSEnd(th,nThreads);
	for (vertexPtr_LGS_iterator it=vertexLGSBegin(th,nThreads);it!=itv_end;++it)
	  {
	    const Coord *c = it->getCoordsConstPtr();
	    for (int i=0;i<NDIM;++i) 
	      {
		Coord checked=geometry->checkCoordConsistency(c[i],ref[i],i);
		if (min[i]>checked) min[i]=checked;
		if (max[i]<checked) max[i]=checked;
	      }	    
	  }

#pragma omp critical
	{
	  for (int i=0;i<NDIM;++i) 
	    {
	      if (xMin[i]>min[i]) xMin[i]=min[i];
	      if (xMax[i]<max[i]) xMax[i]=max[i];
	    }	  
	}
      }
  }

  template <class S, int D=S::NDIM_W, class TT=Coord>
  TT computeVolume(const S *s, bool withSign=false) const
  {
    TT base[S::NVERT-1][D];
    getBaseVectors<S,TT,D>(s,base);   
    if (withSign)
      return SimplexVolumeT<S::NVERT-1,D,TT>::compute_S(base);
    else
      return SimplexVolumeT<S::NVERT-1,D,TT>::compute(base);
  }

  template <class S, class TT=Coord>
  TT computeProjectedVolume(const S *s, bool withSign=false) const
  {
    return computeVolume<S,S::NDIM,TT>(s,withSign);  
  }

  template <class S, int D=S::NDIM_W, class TT=Coord>
  TT computeVolume2(const S *s, bool withSign=false) const
  {
    Coord base[S::NVERT-1][D];
    getBaseVectors<S,TT,D>(s,base);    
    if (withSign)
      return SimplexVolumeT<S::NVERT-1,D,TT>::compute_2S(base);
    else
      return SimplexVolumeT<S::NVERT-1,D,TT>::compute_2(base);
  }
 
  template <class S, class TT=Coord>
  TT computeProjectedVolume2(const S *s, bool withSign=false) const
  {
    return computeVolume2<S,S::NDIM,TT>(s,withSign);  
  }

  // NOTE: the reference vertex for PBC is vertex 0 
  // SH is a pointer to a simplex or a facet handle
  template <class SH, int D>
  void computeBoundingBox(const SH s, Coord bBox[2][D]) const
  {
    const Coord *refCoords = s->getVertex(0)->getCoordsConstPtr();

    for (int j=0;j<D;++j) 
      bBox[0][j]=bBox[1][j]=refCoords[j];

    for (int j=1;j<s->NVERT;++j)
      {
	Vertex *v=s->getVertex(j);
	const Coord *coords=v->getCoordsConstPtr();
	for (int k=0;k<D;++k)
	  {
	    // We have to be a bit carefull for periodic boundaries !
	    double c = geometry->checkCoordConsistency(coords[k],refCoords[k],k);
	    
	    if (bBox[0][k]>c)
	      bBox[0][k]=c;
	    if (bBox[1][k]<c)
	      bBox[1][k]=c;
	  }	
      }
  }
  
  template <class S, typename CT, typename IT, int D=S::NDIM_W>
  long generateUniformSample(const S *s, CT N, IT coordsOut) const
  {
    Coord base[S::NVERT-1][D];
    getBaseVectors<S,Coord,D>(s,base);  
    const Coord *p0=s->getVertex(0)->getCoordsConstPtr();

    return SimplexUniformSamplerT<S::NVERT-1,D>::
      generateRegular(p0,base,N,coordsOut);    
  }

  template <class S, class V, typename CT, typename IT1, typename IT2, int D=S::NDIM_W>
  long generateUniformSample(const S *s, V *vertexValues, 
			     CT N, IT1 coordsOut, IT2 valueOut) const
  {
    Coord base[S::NVERT-1][D];
    getBaseVectors<S,Coord,D>(s,base);  
    const Coord *p0=s->getVertex(0)->getCoordsConstPtr();
   
    return SimplexUniformSamplerT<S::NVERT-1,D>::
      generateRegular(p0,base,vertexValues,N,coordsOut,valueOut);
  }

  template <class S, typename CT, typename IT, int D=S::NDIM_W>
  long generateRandomSample(const S *s, CT N, IT coordsOut, 
			    typename DRand48_rWapper::RandData *randSeed) const
  {
    Coord base[S::NVERT-1][D];
    getBaseVectors<S,Coord,D>(s,base);  
    const Coord *p0=s->getVertex(0)->getCoordsConstPtr();

    return SimplexUniformSamplerT<S::NVERT-1,D>::
      generateRandom(p0,base,N,coordsOut,randSeed);    
  }

  template <class S, class V, typename CT, typename IT1, typename IT2, int D=S::NDIM_W>
  long generateRandomSample(const S *s, V *vertexValues, 
			    CT N, IT1 coordsOut, IT2 valueOut, 
			    typename DRand48_rWapper::RandData *randSeed) const
  {
    Coord base[S::NVERT-1][D];
    getBaseVectors<S,Coord,D>(s,base);  
    const Coord *p0=s->getVertex(0)->getCoordsConstPtr();
   
    return SimplexUniformSamplerT<S::NVERT-1,D>::
      generateRandom(p0,base,vertexValues,N,coordsOut,valueOut,randSeed);
  }


  const MpiCommunication *getMpiCom() const
  {
    return mpiCom;
  }
  
  // compare two segment handles in a consistent way, even accros different processes
  bool compareSegmentHandlesLess(const SegmentHandle &seg1, const SegmentHandle &seg2, bool useSimplex=false) const
  {
    double l1=geometry->template 
      distance2<Coord,NDIM_W>(seg1->getVertex(0)->getCoordsConstPtr(),
			      seg1->getVertex(1)->getCoordsConstPtr());
    double l2=geometry->template 
      distance2<Coord,NDIM_W>(seg2->getVertex(0)->getCoordsConstPtr(),
			      seg2->getVertex(1)->getCoordsConstPtr());

    if (l1<l2) return true;
    else if (l2<l1) return false;

    if (useSimplex)
      {
	Simplex *s1 = seg1->getSimplex();
	Simplex *s2 = seg2->getSimplex();
	if ((s1!=NULL)&&(s2!=NULL))
	  {
	    l1 = computeVolume(s1);
	    l2 = computeVolume(s2);
	    
	    if (l1<l2) return true;
	    else if (l2<l1) return false;
	    
	    GlobalIdentity g1 = s1->getGlobalIdentity(mpiCom->rank());
	    GlobalIdentity g2 = s2->getGlobalIdentity(mpiCom->rank());
	    if (g1<g2) return true;
	    else if (g2<g1) return false;
	  }
      }

    GlobalIdentity g1a = seg1->getVertex(0)->getGlobalIdentity();
    GlobalIdentity g1b = seg1->getVertex(1)->getGlobalIdentity();
    if (g1b<g1a) std::swap(g1a,g1b);
    GlobalIdentity g2a = seg2->getVertex(0)->getGlobalIdentity();
    GlobalIdentity g2b = seg2->getVertex(1)->getGlobalIdentity();
    if (g2b<g2a) std::swap(g2a,g2b);
    
    if (g1a<g2a) return true;
    else if (g2a<g1a) return false;
 
    if (g1b<g2b) return true;
    else if (g2b<g1b) return false;

    return false;
  }

  std::vector<ShadowVertex*> getShadowVerticesArray()
  {
    std::vector<ShadowVertex*> result(shadowVertexPool.getUsedCount());
    if (result.size()==0) return result;
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const shadowVertexPtr_iterator it_end=shadowVertexEnd();
	for (shadowVertexPtr_iterator it=shadowVertexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }
    

  std::vector<ShadowSimplex*> getShadowSimplicesArray()
  {
    std::vector<ShadowSimplex*> result(shadowSimplexPool.getUsedCount());
    if (result.size()==0) return result;
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const shadowSimplexPtr_iterator it_end=shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=shadowSimplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }

  std::vector<GhostVertex*> getGhostVerticesArray()
  {
    std::vector<GhostVertex*> result(ghostVertexPool.getUsedCount());
    if (result.size()==0) return result;
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const ghostVertexPtr_iterator it_end=ghostVertexEnd();
	for (ghostVertexPtr_iterator it=ghostVertexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }

  std::vector<GhostSimplex*> getGhostSimplicesArray()
  {
    std::vector<GhostSimplex*> result(ghostSimplexPool.getUsedCount());
    if (result.size()==0) return result;
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const ghostSimplexPtr_iterator it_end=ghostSimplexEnd();
	for (ghostSimplexPtr_iterator it=ghostSimplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }
  
  void resetSimplicesCache(int nThreads=glb::num_omp_threads)
  {
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<nThreads;i++)
      {	
	const simplexPtr_iterator it_end=simplexEnd();
	for (simplexPtr_iterator it=simplexBegin(i,nThreads);it!=it_end;++it)
	  it->resetCache();
      }  
  }

  void resetGhostSimplicesCache(int nThreads=glb::num_omp_threads)
  {
    if (getNGhostSimplices()==0) return;
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<nThreads;i++)
      {	
	const ghostSimplexPtr_iterator it_end=ghostSimplexEnd();
	for (ghostSimplexPtr_iterator it=ghostSimplexBegin(i,nThreads);it!=it_end;++it)
	  {
	    memset(&it->cache,0,sizeof(it->cache));
	  }
      }      
  }

  void setSimplicesCacheL(long val)
  {
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {	
	const simplexPtr_iterator it_end=simplexEnd();
	for (simplexPtr_iterator it=simplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    it->cache.l=val;	   
	  }
      }      
  }

protected: 
  MpiCommunication *mpiCom;      

  SimplexPool simplexPool;  
  VertexPool vertexPool;
  GhostSimplexPool ghostSimplexPool;
  GhostVertexPool ghostVertexPool;
  ShadowSimplexPool shadowSimplexPool; 
  ShadowVertexPool shadowVertexPool;
 
  std::vector<double> x0;
  std::vector<double> delta;
  
  CellDataFunctors cellDataFunctors;

  GeometricProperties *geometry;
  //MockSimplexData *volumeFunctor;
  //MockSimplexData *projectedVolumeFunctor;
  bool initialized;
  /*
  template <class M>
  class ProjectedVolumeT
  {
  public:
    ProjectedVolumeT(const M* m):mesh(m){}
    double get(const Simplex *s) const
    {
      return mesh->computeProjectedVolume(s);
    }
  private:
    const M* mesh;
  };

  template <class M>
  class VolumeT
  {
  public:
    VolumeT(const M* m):mesh(m){}
    double get(const Simplex *s) const
    {
      return mesh->computeVolume(s);
    }
  private:
    const M* mesh;   
  };

  VolumeT<MyType> *volume;
  ProjectedVolumeT<MyType> *projectedVolume;
  */

  template <class W>
  void write(W *writer)
  {    
    writer->writeHeader(classHeader(),classVersion());  
      
    vertexPool.serialize(writer);
    ghostVertexPool.serialize(writer);
    shadowVertexPool.serialize(writer);
	    
    simplexPool.serialize(writer);
    ghostSimplexPool.serialize(writer);
    shadowSimplexPool.serialize(writer);

    Tree::write(writer);
  }

  void resetVerticesTags(int nThreads=glb::num_omp_threads)
  {
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<nThreads;i++)
      {
	const vertexPtr_iterator it_end=vertexEnd();
	for (vertexPtr_iterator it=vertexBegin(i,nThreads);it!=it_end;++it)
	  {
	    it->cleanTags();
	  }
      }
  }

  // This is multithreaded, but it is the caller's responsibility to check
  // that there are no refinement conflicts.
  // sharedVertices will contain vertices that are shared and
  // whose globalIdentity should therefore be updated (caller's responsibility).
  // To facilitate this process, the incomplete shared vertices are given the globalIdentity of 
  // the remote simplex that created them (i.e. refinedSegment->getSimplex())
  // Note that splitSegments only refine local simplices, not ghosts/shadows
  // NOTE: this will set simplices cache ptr for split simplices to point to their partner
  template <class C, class M, typename TT, class S>
  long splitSegments(const C &refinedSegments,
		     M &newSharedVerticesMap,
		     std::vector<Vertex*> &newVertices,
		     std::vector<Simplex*> &newSimplices,
		     std::vector<TT> &nSimplicesCum,
		     S *solver, int nThreads)
		     
  {    
    // resetSimplicesCache(); Not always necessary, so caller may want to do that ...
    long nRefined=refinedSegments.size();
    if (nRefined==0) return 0;    
 
    // count how many simplices will be refined for each segment
    //std::vector<unsigned long> nSimplicesCum(refinedSegments.size()+1);
    nSimplicesCum.assign(refinedSegments.size()+1,0);
  
    std::vector<unsigned char> segmentIsShared(refinedSegments.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<refinedSegments.size();i++)
      segmentIsShared[i]=
	refinedSegments[i].second->countLocalSimplices(nSimplicesCum[i+1]); 
 
    /*
    for (unsigned long i=0;i<refinedSegments.size();i++)
      if ((refinedSegments[i].second->getVertex(0)->getGlobalIdentity()==GlobalIdentity(0,19))||
	  (refinedSegments[i].second->getVertex(1)->getGlobalIdentity()==GlobalIdentity(0,375)))
	{
	  refinedSegments[i].second->getSimplex()->template print<LOG_STD_ALL>("TOREFINE");
	  refinedSegments[i].second->getSimplex()->template print<LOG_STD_ALL>("TOREFINE");
	  refinedSegments[i].second->template print<LOG_STD_ALL>();
	}
    */

    nSimplicesCum[0]=0;
    for (unsigned long i=0;i<refinedSegments.size();i++)
      nSimplicesCum[i+1]+=nSimplicesCum[i];

    // allocate new elements
    // std::vector<Vertex*> newVertices;    
    std::vector<Segment*> newSegments;
    std::vector<Node*> newNodes;
    //std::vector<Simplex*> newSimplices;   
    
    newVertices.reserve(refinedSegments.size());
    newSimplices.reserve(nSimplicesCum.back());
    newNodes.reserve(nSimplicesCum.back());
      
    vertexPool.reserve(refinedSegments.size());
    simplexPool.reserve(nSimplicesCum.back());
    Tree::nodePool.reserve(nSimplicesCum.back());    

    std::vector<Coord> dummyCoord(NDIM_W,0);
    for (unsigned long i=0;i<refinedSegments.size();i++)
      {	  
	// this happens when all the refined simplices are ghosts and/or shadows
	if (nSimplicesCum[i] == nSimplicesCum[i+1]) 
	  {
	    // newVertices are adressed by indexing newVertices
	    // array so we must still insert a dummy vertex
	    // This is also used to spot segments
	    // that do not refine anything locally
	    newVertices.push_back(NULL);
	    continue;
	  }

	Vertex *v;
	Simplex *s;
	Node *n;
	Simplex *rs=refinedSegments[i].second->getSimplex();
	
	vertexPool.pop(&v);
	v->set(&dummyCoord[0],getNCells(0)-1);
	v->setGlobalIdentity(mpiCom->rank(),v->getLocalIndex());
	v->setGeneration(refinedSegments[i].second->getGeneration(),0);
	if (segmentIsShared[i]) 
	  {
	    v->setSharedF();
	    // the key is the (correct) globalIdentity of the simplex that caused splitting
	    newSharedVerticesMap.insert(std::make_pair(rs->getGlobalIdentity(mpiCom->rank()).get(),v));
	  }

	newVertices.push_back(v);

	for (TT j=nSimplicesCum[i];j<nSimplicesCum[i+1];j++)
	  {
	    simplexPool.pop(&s);
	    s->setLocalIndex(getNCells(NDIM)-1);
	    newSimplices.push_back(s);
	    Tree::nodePool.pop(&n);
	    newNodes.push_back(n);	    
	  }	
      }

    nRefined = refinedSegments.size();

    const long rBufferSize=Vertex::template 
      getSimplexRefineBufferSize<MyType,SegmentHandle>();

    std::vector<char> simplexBuffer;
    simplexBuffer.resize(newSimplices.size()*rBufferSize);
   
    // FIXME : We could go faster in what follows by storing the simplices 
    // in 'lst' computed in the first loop and reuse them in the second loop.
    // But this is not a critical part of the code so we don t really care ...

    // Now that we are sure there is no conflict left, do the refinement.  
    // First we compute the new vertices coordinates, generation and data
#pragma omp parallel for num_threads(nThreads) schedule(dynamic,1)
    for (unsigned long i=0;i<refinedSegments.size();i++)
      {
	Vertex *v=newVertices[i];
	if (v==NULL) continue;
	const SegmentHandle &h=refinedSegments[i].second;

	int nSimplices=nSimplicesCum[i+1]-nSimplicesCum[i];
	Simplex *lst[nSimplices];
	char* buffer=&simplexBuffer[nSimplicesCum[i]*rBufferSize];
	//Simplex **s= &newSimplices[nSimplicesCum[i]];

	//v->setGeneration(h->getGeneration(),0);
	
	segment_circulator ci_end=h->getCirculator();
	segment_circulator ci=ci_end;
	int j=0;
	do
	  {
	    if (ci->isLocal())
	      lst[j++]=*ci;
	  } while ((++ci)!=ci_end);

	Vertex::Data::template onRefineVertex<MyType,SegmentHandle,Vertex,Simplex>
	  (this,h,v,lst,nSimplices,buffer);
	//solver->onRefineVertex(h,v,&lst[0],nSimplices);
      }
    
    // then split the simplices and recompute the data
#pragma omp parallel for num_threads(nThreads) schedule(dynamic,1)
    for (unsigned long i=0;i<refinedSegments.size();i++)
      {
	Vertex *v=newVertices[i];

	// this happens when all the refined simplices are ghosts and/or shadows
	// In that case, there is nothing to do locally ...
	if (v==NULL) continue;
	
	const SegmentHandle &h=refinedSegments[i].second;		  
	Node **n = &newNodes[nSimplicesCum[i]];
	Simplex **s= &newSimplices[nSimplicesCum[i]];
	char* buffer=&simplexBuffer[nSimplicesCum[i]*rBufferSize];
	int nSimplices=nSimplicesCum[i+1]-nSimplicesCum[i];
	Simplex *lst[nSimplices];	  	  
	
	// v->setGeneration(h->getGeneration(),0);		
	// solver->onRefineVertex(h,v,lst,nSimplices);

	// We cannot split simplices within the loop without messing-up
	// the segment circulator
	segment_circulator ci_end=h->getCirculator();
	segment_circulator ci=ci_end;	
	int j=0;
	do
	  {
	    if (ci->isLocal())
	      lst[j++]=*ci;
	  } while ((++ci)!=ci_end);
	
	// split every simplex that contains segment h=refinedSegments[i]
	for (j=0;j<nSimplices;j++)
	  lst[j]->splitSegment(h,n[j],v,s[j]);	    
	 
	Simplex::Data::template onRefineSimplices<MyType,SegmentHandle,Vertex,Simplex>
	  (this,h,v,lst,s,nSimplices,buffer);
	//solver->onRefineSimplices(h,v,lst,s,nSimplices);
      }
    
    // FIXME: getRoot is ~log(N) where N is the number of simplices ...
    // We do it here because it is not thread safe !
    for (TT i=0;i<nSimplicesCum.back();i++)
      newSimplices[i]->getRoot()->increaseWeight();

    return newVertices.size();
  }

  // Split ONE segment, this is not thread safe ...
  // + UNTESTED !!!
  /*
  void splitSegment(SegmentHandle &h)
  {
    std::vector<Simplex*> lst;
    int j=0;
    Vertex *v;
    vertexPool.pop(v);  
    
    segment_circulator ci_end=h->getSimplex()->getSegmentCirculator(h);
    segment_circulator ci=ci_end;	  
    do
      {
	lst[j++]=*ci;
      } while ((++ci)!=ci_end);

    Node* nodes[2];
    Simplex* newSimplices[2*lst.size()];

    for (j=0;j<lst.size();j++)
      {
	Tree::nodesPool.pop(&nodes[0]);
	Tree::nodesPool.pop(&nodes[1]);
	simplexPool.pop(&newSimplices[2*j+0]);
	simplexPool.pop(&newSimplices[2*j+1]);
	lst[j]->splitSegment(h,nodes,v,&newSimplices[2*j]);
      }

    // fix neighbors
    for (int i=0;i<lst.size()*2;i+=2)
      {
	newSimplices[i]->fixNeighborsAfterSplitting();
	Simplex *other=static_cast<Simplex*>(newSimplices[i]->getOther());
	other->fixNeighborsAfterSplitting();
      }

    // fix neighbors's neighbors
    Simplex **start = &newSimplices[0];
    Simplex **stop = &newSimplices[2*lst.size()];
    for (Simplex **s=start;s!=stop;++s)
      {
	int id=(*s)->getVertexIndex(v);
	Simplex *nei=(*s)->getNeighbor(id);      
	if (nei!=NULL)
	  {
	    if (!nei->isShadow()) 
	      {
		id=nei->getOppositeVertexIndex(*s);
		if (id<0) exit(0);
		nei->neighbors[id]=*s;		  
	      }
	  }
	Simplex *other=static_cast<Simplex *>(*s)->getOther();
	id=other->getVertexIndex(v);
	nei=other->getNeighbor(id);
	if (nei!=NULL)
	  {
	    if (!nei->isShadow())
	      {
		id=nei->getOppositeVertexIndex(other);
		if (id<0) exit(0);
		nei->neighbors[id]=other;
	      }
	  }
      }
  }
  */

  
  // Used to correct the simplices neighbors that have became wrong after splitting.
  // The newly created simplices associated to new vertex newVertices[i] are 
  // thoses in the range (newSimplices[nSimplicesCum[i]],newSimplices[nSimplicesCum[i+1]](
  // When skipShared is true, shadows and ghosts are ignored
  template <class TT>
  void fixSimplicesNeighborsAfterSplitting(std::vector<Vertex*> &newVertices,
					   std::vector<Simplex*> &newSimplices, 
					   std::vector<TT> &nSimplicesCum, 
					   bool skipShared=true)
  {
    if (newVertices.size()==0) return;    

    // fix the neighbors that have been split
    //FIXME #pragma omp parallel for
    for (unsigned long i=0;i<newSimplices.size();i++)
      {
	// get the partner of the new simplex
	Simplex *other=newSimplices[i]->
	  getNeighbor(newSimplices[i]->cache.c[0]); 

	newSimplices[i]->fixNeighborsAfterSplitting(skipShared);	
	other->fixNeighborsAfterSplitting(skipShared);
      }
    
    // as well as the neighbors that haven't been split.
    // Now that split neighbors are fixed, for each split simplex,
    // only the simplices opposite to the new vertex need fixing. 
    //FIXME: omp here, but nei->setNeighbor(id,cur); is critical ...
    // #pragma omp parallel for
    for (unsigned long i=0;i<newVertices.size();i++)
      {
	Simplex **start = &newSimplices[nSimplicesCum[i]];
	Simplex **stop = &newSimplices[nSimplicesCum[i+1]];
	for (Simplex **s=start;s!=stop;++s)
	  {
	    Simplex *cur=*s;
	    Simplex *other=(*s)->getNeighbor((*s)->cache.c[0]); //get the partner of cur

	    int id=cur->getVertexIndex(newVertices[i]);
	    Simplex *nei=cur->getNeighbor(id);      
	    if (nei!=NULL)
	      {
		bool isShared = nei->isShadowOrGhost();
		if (!skipShared)
		  {
		    id=nei->getOppositeVertexIndex(cur);
		    if (id<0) 
		      {
			PRINT_SRC_INFO(LOG_ERROR);
			glb::console->print<LOG_ERROR>("Disconnected split simplex neighbor.\n");		
			glb::console->print<LOG_ERROR>("That's embarrassing (1A)\n");
			cur->template print<LOG_ERROR>("CUR:");
			nei->template print<LOG_ERROR>("NEI:");
			newVertices[i]->template print<LOG_ERROR>();	
			exit(-1);
		      }
		    // if (nei->getGlobalIdentity(0).id() == 2816)
		    //   {
			
		    //   	glb::console->print<LOG_WARNING>("\nERROR ON CULP 1\n");
		    //  	nei->template print<LOG_WARNING>(); 
		    //   	cur->template print<LOG_WARNING>();

		    //   }
		    nei->setNeighbor(id,cur);
		    //nei->neighbors[id]=cur;	
		  }
		else if (!isShared)
		  {
		    id=nei->getOppositeVertexIndex(cur);
		    if (id<0) 
		      {
			PRINT_SRC_INFO(LOG_ERROR);
			glb::console->print<LOG_ERROR>("Disconnected split simplex neighbor!\n");
			glb::console->print<LOG_ERROR>("That's embarrassing (1B)\n");
			exit(-1);
		      }
		    nei->setNeighbor(id,cur);
		    //nei->neighbors[id]=cur;		  
		  }
	      }
	    
	    id=other->getVertexIndex(newVertices[i]);
	    nei=other->getNeighbor(id);
	    if (nei!=NULL)
	      {
		bool isShared = nei->isShadowOrGhost(); 
		if (!skipShared)
		  {
		    id=nei->getOppositeVertexIndex(other);
		    if (id<0)
		      {
			PRINT_SRC_INFO(LOG_ERROR);
			glb::console->print<LOG_ERROR>("Disconnected split simplex neighbor!\n");		
			glb::console->print<LOG_ERROR>("That's embarrassing (2A)\n");
			cur->template print<LOG_ERROR>("CUR:");
			nei->template print<LOG_ERROR>("NEI:");
			newVertices[i]->template print<LOG_ERROR>();	
			exit(-1);
		      }
		    // if (nei->getGlobalIdentity(0).id() == 2816)
		    //     {
			
		    //  	glb::console->print<LOG_WARNING>("\nERROR ON CULP 2\n");
		    //   	nei->template print<LOG_WARNING>(); 
		    // 	nei->getNeighbor(0)->template print<LOG_WARNING>(); 
		    //   	other->template print<LOG_WARNING>();
		    //    }
		    nei->setNeighbor(id,other);
		    //nei->neighbors[id]=other;
		  }
		else if (!isShared)
		  {
		    id=nei->getOppositeVertexIndex(other);
		    if (id<0)
		      {
			PRINT_SRC_INFO(LOG_ERROR);
			glb::console->print<LOG_ERROR>("Disconnected split simplex neighbor!\n");
			glb::console->print<LOG_ERROR>("That's embarrassing (2B)\n");
			exit(-1);
		      }
		    nei->setNeighbor(id,other);
		    //nei->neighbors[id]=other;
		  }
	      }
	  }
      }
  }
  
  template <class LOG>
  bool checkConsistency(bool checkDuplicates=false)
  {
    bool success = true;
    const simplexPtr_iterator its_end=simplexEnd();
    for (simplexPtr_iterator it=simplexBegin();it!=its_end;++it)
      {
	Simplex *s=*it;
	if (!s->template checkConsistency<LOG>(false))
	  {
	    glb::console->print<LOG>("ERROR checking simplex !\n");	    
	    s->template checkConsistency<LOG>();
	    s->template print<LOG>();
	    success=false;
	  }
	for (int i=0;i<Simplex::NNEI;++i)
	  {
	    Simplex *ns=(*it)->getNeighbor(i);
	    if (ns==NULL) continue;
	    bool res;
	    if (ns->isShadow())
	      res=shadowSimplexPool.isFree(ns->getAsShadow());
	    else if (ns->isGhost())
	      res=ghostSimplexPool.isFree(ns->getAsGhost());
	    else
	      res=simplexPool.isFree(ns);
	    if (res==true)
	      {
		glb::console->print<LOG>("ERROR checking ghost simplex !\n");
		glb::console->print<LOG>("Neighbor @%d was freed !\n",i);
		s->template print<LOG>();
		success=false;
	      }
	  }	
      }

    const ghostSimplexPtr_iterator itgs_end=ghostSimplexEnd();
    for (ghostSimplexPtr_iterator it=ghostSimplexBegin();it!=itgs_end;++it)
      {
	GhostSimplex *s=*it;
	if (!s->template checkConsistency<LOG>(false))
	  {
	    glb::console->print<LOG>("ERROR checking ghost simplex !\n");	    
	    s->template checkConsistency<LOG>();
	    s->template print<LOG>();
	    success=false;
	  }
	for (int i=0;i<Simplex::NNEI;++i)
	  {
	    Simplex *ns=(*it)->getNeighbor(i);
	    if (ns==NULL) continue;
	    bool res;
	    if (ns->isShadow())
	      res=shadowSimplexPool.isFree(ns->getAsShadow());
	    else if (ns->isGhost())
	      res=ghostSimplexPool.isFree(ns->getAsGhost());
	    else
	      res=simplexPool.isFree(ns);
	    if (res==true)
	      {
		glb::console->print<LOG>("ERROR checking ghost simplex !\n");
		glb::console->print<LOG>("Neighbor @%d was freed !\n",i);
		s->template print<LOG>();
		success=false;
	      }
	  }
      }

    const shadowSimplexPtr_iterator itss_end=shadowSimplexEnd();
    for (shadowSimplexPtr_iterator it=shadowSimplexBegin();it!=itss_end;++it)
      {
	ShadowSimplex *s=*it;
	if (!s->template checkConsistency<LOG>(false))
	  {
	    glb::console->print<LOG>("ERROR checking shadow simplex !\n");	    
	    s->template checkConsistency<LOG>();
	    s->template print<LOG>();
	    success=false;
	  }
	for (int i=0;i<Simplex::NNEI;++i)
	  {
	    Simplex *ns=(*it)->getNeighbor(i);
	    if (ns==NULL) continue;
	    bool res;
	    if (ns->isShadow())
	      res=shadowSimplexPool.isFree(ns->getAsShadow());
	    else if (ns->isGhost())
	      res=ghostSimplexPool.isFree(ns->getAsGhost());
	    else
	      res=simplexPool.isFree(ns);
	    if (res==true)
	      {
		glb::console->print<LOG>("ERROR checking shadow simplex !\n");
		glb::console->print<LOG>("Neighbor @%d was freed !\n",i);
		s->template print<LOG>();
		success=false;
	      }
	  }
      }

    if (checkDuplicates)
      {
	const long myRank=mpiCom->rank();

	typedef typename my_dense_hash<GlobalIdentityValue,Simplex*>::type 
	  GidSimplexHash;
	typedef typename GidSimplexHash::iterator GidSimplexHash_it;

	typedef typename my_dense_hash<GlobalIdentityValue,Vertex*>::type 
	  GidVertexHash;
	typedef typename GidVertexHash::iterator GidVertexHash_it;

	GidSimplexHash simplexHash;
	GidVertexHash vertexHash;
	set_hash_empty_key(simplexHash,GlobalIdentity::empty.get());
	set_hash_empty_key(vertexHash,GlobalIdentity::empty.get());
	simplexHash.rehash(getNCellsTotal(NDIM)/
			   (simplexHash.max_load_factor()-0.01));
	vertexHash.rehash(getNCellsTotal(0)/
			  (vertexHash.max_load_factor()-0.01));

	for (simplexPtr_iterator it=simplexBegin();it!=its_end;++it)
	  {
	    std::pair<GidSimplexHash_it,bool> res = 
	      simplexHash.insert(std::make_pair(it->getGlobalIdentity(myRank).get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate Simplex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }
	for (ghostSimplexPtr_iterator it=ghostSimplexBegin();it!=itgs_end;++it)
	  {
	    std::pair<GidSimplexHash_it,bool> res = 
	      simplexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate Simplex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }
	for (shadowSimplexPtr_iterator it=shadowSimplexBegin();it!=itss_end;++it)
	  {
	    std::pair<GidSimplexHash_it,bool> res = 
	      simplexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate Simplex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }

	const vertexPtr_iterator itv_end=vertexEnd();
	for (vertexPtr_iterator it=vertexBegin();it!=itv_end;++it)
	  {
	    std::pair<GidVertexHash_it,bool> res = 
	      vertexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate Vertex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }
	const ghostVertexPtr_iterator itgv_end=ghostVertexEnd();
	for (ghostVertexPtr_iterator it=ghostVertexBegin();it!=itgv_end;++it)
	  {
	    std::pair<GidVertexHash_it,bool> res = 
	      vertexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate Vertex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }
	const shadowVertexPtr_iterator itsv_end=shadowVertexEnd();
	for (shadowVertexPtr_iterator it=shadowVertexBegin();it!=itsv_end;++it)
	  {
	    std::pair<GidVertexHash_it,bool> res = 
	      vertexHash.insert(std::make_pair(it->getGlobalIdentity().get(),*it));
	    if (!res.second)
	      {
		glb::console->print<LOG>("Duplicate vertex :\n");
		(*it)->template print<LOG>();
		(res.first->second)->template print<LOG>();
		success=false;
	      }
	  }
      }

    return success;
  }

  template <class SimplexFunctor, class VertexFunctor,
	    class SPU, class GSPU, class SSPU,
	    class VPU, class GVPU, class SVPU>
  void sort(const SimplexFunctor &sf, 
	    const VertexFunctor &vf,
	    SPU &lspu, GSPU &gspu, SSPU &sspu,
	    VPU &lvpu, GVPU &gvpu, SVPU &svpu,	   
	    int nThreads=glb::num_omp_threads)
  {        
    glb::console->printFlush<LOG_INFO>("Sorting simplices ... ");   

    lspu=simplexPool.sort(sf,nThreads); 
    gspu=ghostSimplexPool.sort(sf,nThreads); 
    sspu=shadowSimplexPool.sort(sf,nThreads); 

    glb::console->printFlush<LOG_INFO>("done.\n");   

    glb::console->printFlush<LOG_INFO>("Updating simplex pointers ... ");   
#pragma omp parallel for num_threads(nThreads) schedule(static,1)
    for (long i=0;i<nThreads;i++)
      {
	const simplexPtr_LGS_iterator it_end=simplexLGSEnd();
	for (simplexPtr_LGS_iterator it=simplexLGSBegin(i,nThreads);
	     it!=it_end;++it)
	  (*it)->updateSimplicesPointers(lspu,gspu,sspu);	
      } 
    glb::console->printFlush<LOG_INFO>("done.\n");

    // Now that the simplices have been sorted, we need a function to 
    // sort the vertices accordingly !

    glb::console->printFlush<LOG_INFO>("Sorting vertices ... ");   

    lvpu=vertexPool.sort(vf,nThreads); 
    gvpu=ghostVertexPool.sort(vf,nThreads); 
    svpu=shadowVertexPool.sort(vf,nThreads); 

    glb::console->printFlush<LOG_INFO>("done.\n");

    glb::console->printFlush<LOG_INFO>("Updating vertex pointers ... ");   
#pragma omp parallel for num_threads(nThreads) schedule(static,1) 
    for (long i=0;i<nThreads;i++)
      {
	const simplexPtr_LGS_iterator it_end=simplexLGSEnd();
	for (simplexPtr_LGS_iterator it=simplexLGSBegin(i,nThreads);
	     it!=it_end;++it)
	  (*it)->updateVerticesPointers(lvpu,gvpu,svpu);	
      } 
    glb::console->printFlush<LOG_INFO>("done.\n");
  }
  

  /*
  std::vector<ShadowSimplex*> getShadowSimplicesArray(bool sorted=false)
  {
    std::vector<ShadowSimplex*> result;

    const long nTot=shadowSimplexPool.getUsedCount();   
    std::list<ShadowSimplex*> remainder;
    std::list<unsigned long> freeId;
    typedef typename std::list<ShadowSimplex*>::iterator remainderIt;
    typedef typename std::list<unsigned long>::iterator freeIdIt;

    result.resize(nTot);

#pragma omp parallel for
    for (unsigned long i=0;i<glb::num_omp_threads;i++)
      {
	std::list<ShadowSimplex*> localRemainder;
	unsigned long j=i;
	const shadowSimplexPtr_iterator it_end=shadowSimplexEnd();
	for (shadowSimplexPtr_iterator it=shadowSimplexBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {
	    if (j>=nTot)
	      localRemainder.push_back(*it);
	    else
	      result[j]=(*it);	
	    
	    j+=glb::num_omp_threads;	      
	  }

	// this is for when the number of vertices on each proc is not balanced ...
	if (j>nTot+i)
	  {
#pragma omp critical
	    remainder.splice(remainder.end(),localRemainder,localRemainder.begin(),localRemainder.end());
	  }
	else if (j<nTot)
	  {
#pragma omp critical
	    freeId.push_back(j);
	  }
      }
    
    remainderIt it = remainder.begin();
    for (freeIdIt free_it=freeId.begin(); free_it !=freeId.end(); ++free_it)
      {
	for (unsigned int j=(*free_it); j<nTot; j+=glb::num_omp_threads)
	  {
	    result[j]=(*it);
	    ++it;
	  }
      }    
    
    if (sorted)
      {
#ifdef USE_GNU_PSORT
	__gnu_parallel::sort(result.begin(), result.end(), hlp::comparePointersToLess<ShadowSimplex>());
#else
	std::sort(result.begin(),result.end(), hlp::comparePointersToLess<ShadowSimplex>());
#endif
      }
    
    return result;
  }
  */


  // iterators definitions
public:  

  typedef typename SimplexPool::iterator simplexPtr_iterator;
  typedef typename VertexPool::iterator vertexPtr_iterator;  

  typedef typename SimplexPool::const_iterator const_simplexPtr_iterator;
  typedef typename VertexPool::const_iterator const_vertexPtr_iterator;

  typedef typename GhostSimplexPool::iterator ghostSimplexPtr_iterator;
  typedef typename GhostVertexPool::iterator ghostVertexPtr_iterator;
  typedef typename GhostSimplexPool::const_iterator const_ghostSimplexPtr_iterator;
  typedef typename GhostVertexPool::const_iterator const_ghostVertexPtr_iterator;

  typedef typename ShadowSimplexPool::iterator shadowSimplexPtr_iterator;
  typedef typename ShadowVertexPool::iterator shadowVertexPtr_iterator;
  typedef typename ShadowSimplexPool::const_iterator const_shadowSimplexPtr_iterator;
  typedef typename ShadowVertexPool::const_iterator const_shadowVertexPtr_iterator;
  
  typedef UnionIterator3T<simplexPtr_iterator,
			  ghostSimplexPtr_iterator,
			  shadowSimplexPtr_iterator> simplexPtr_LGS_iterator;
  typedef UnionIterator3T<const_simplexPtr_iterator,
			  const_ghostSimplexPtr_iterator,
			  const_shadowSimplexPtr_iterator> const_simplexPtr_LGS_iterator;
  typedef UnionIterator3T<vertexPtr_iterator,
			  ghostVertexPtr_iterator,
			  shadowVertexPtr_iterator> vertexPtr_LGS_iterator;
  typedef UnionIterator3T<const_vertexPtr_iterator,
			  const_ghostVertexPtr_iterator,
			  const_shadowVertexPtr_iterator> const_vertexPtr_LGS_iterator;

  typedef UnionIterator2T<simplexPtr_iterator,
			  ghostSimplexPtr_iterator> simplexPtr_LG_iterator;
  typedef UnionIterator2T<const_simplexPtr_iterator,
			  const_ghostSimplexPtr_iterator> const_simplexPtr_LG_iterator;
  typedef UnionIterator2T<vertexPtr_iterator,
			  ghostVertexPtr_iterator> vertexPtr_LG_iterator;
  typedef UnionIterator2T<const_vertexPtr_iterator,
			  const_ghostVertexPtr_iterator> const_vertexPtr_LG_iterator;
  
  simplexPtr_iterator simplexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return simplexPool.begin();
    else
      return simplexPool.begin(delta,stride);
  }
  
  simplexPtr_iterator simplexEnd(int delta=0, int stride=1)
  {
    return simplexPool.end();
  }
  /*
  segmentPtr_iterator segmentBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return segmentPool.begin();
    else
      return segmentPool.begin(delta,stride);
  }
  
  segmentPtr_iterator segmentEnd()
  {
    return segmentPool.end();
  }
  */
  vertexPtr_iterator vertexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return vertexPool.begin();
    else
      return vertexPool.begin(delta,stride);
  }
  
  vertexPtr_iterator vertexEnd(int delta=0, int stride=1)
  {
    return vertexPool.end();
  }

  shadowSimplexPtr_iterator shadowSimplexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return shadowSimplexPool.begin();
    else
      return shadowSimplexPool.begin(delta,stride);
  }
  
  shadowSimplexPtr_iterator shadowSimplexEnd(int delta=0, int stride=1)
  {
    return shadowSimplexPool.end();
  }

  shadowVertexPtr_iterator shadowVertexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return shadowVertexPool.begin();
    else
      return shadowVertexPool.begin(delta,stride);
  }
  
  shadowVertexPtr_iterator shadowVertexEnd(int delta=0, int stride=1)
  {
    return shadowVertexPool.end();
  }

  ghostVertexPtr_iterator ghostVertexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return ghostVertexPool.begin();
    else
      return ghostVertexPool.begin(delta,stride);
  }
  
  ghostVertexPtr_iterator ghostVertexEnd(int delta=0, int stride=1)
  {
    return ghostVertexPool.end();
  }

  ghostSimplexPtr_iterator ghostSimplexBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return ghostSimplexPool.begin();
    else
      return ghostSimplexPool.begin(delta,stride);
  }
  
  ghostSimplexPtr_iterator ghostSimplexEnd(int delta=0, int stride=1)
  {
    return ghostSimplexPool.end();
  }

  vertexPtr_LGS_iterator vertexLGSBegin(int delta=0, int stride=1)
  {
    long d1 = delta - (getNVertices()%stride);
    if (d1<0) d1+=stride;
    long d2 = d1 - (getNGhostVertices()%stride);
    if (d2<0) d2+=stride;
    //glb::console->print<LOG_DEBUG>("(delta=%d,stride=%d) -> (d1=%ld) (d2=%ld)\n",delta,stride,d1,d2);
    return vertexPtr_LGS_iterator(vertexBegin(delta,stride),
				  vertexEnd(delta,stride),
				  ghostVertexBegin(d1,stride),
				  ghostVertexEnd(d1,stride),
				  shadowVertexBegin(d2,stride),
				  shadowVertexEnd(d2,stride));
  }

  vertexPtr_LGS_iterator vertexLGSEnd(int delta=0, int stride=1)
  {
    return vertexPtr_LGS_iterator(vertexEnd(delta,stride),
				  ghostVertexEnd(delta,stride),
				  shadowVertexEnd(delta,stride));
  }

  simplexPtr_LGS_iterator simplexLGSBegin(int delta=0, int stride=1)
  {
    long d1 = delta - (getNSimplices()%stride);
    if (d1<0) d1+=stride;
    long d2 = d1 - (getNGhostSimplices()%stride);
    if (d2<0) d2+=stride;
    return simplexPtr_LGS_iterator(simplexBegin(delta,stride),
				   simplexEnd(delta,stride),
				   ghostSimplexBegin(d1,stride),
				   ghostSimplexEnd(d1,stride),
				   shadowSimplexBegin(d2,stride),
				   shadowSimplexEnd(d2,stride));
  }

  simplexPtr_LGS_iterator simplexLGSEnd(int delta=0, int stride=1)
  {
    return simplexPtr_LGS_iterator(simplexEnd(delta,stride),
				   ghostSimplexEnd(delta,stride),
				   shadowSimplexEnd(delta,stride));
  }


  vertexPtr_LG_iterator vertexLGBegin(int delta=0, int stride=1)
  {
    long d1 = delta - (getNVertices()%stride);
    if (d1<0) d1+=stride;  
    return vertexPtr_LG_iterator(vertexBegin(delta,stride),vertexEnd(delta,stride),
				 ghostVertexBegin(d1,stride),ghostVertexEnd(d1,stride));
  }

  vertexPtr_LG_iterator vertexLGEnd(int delta=0, int stride=1)
  {
    return vertexPtr_LG_iterator(vertexEnd(delta,stride),ghostVertexEnd(delta,stride));
  }

  simplexPtr_LG_iterator simplexLGBegin(int delta=0, int stride=1)
  {
    long d1 = delta - (getNSimplices()%stride);
    if (d1<0) d1+=stride;    
    return simplexPtr_LG_iterator(simplexBegin(delta,stride),
				  simplexEnd(delta,stride),
				  ghostSimplexBegin(d1,stride),
				  ghostSimplexEnd(d1,stride));
  }

  simplexPtr_LG_iterator simplexLGEnd()
  {
    return simplexPtr_LG_iterator(simplexEnd(),ghostSimplexEnd());
  }


  // CONST iterators

  const_simplexPtr_iterator simplexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return simplexPool.begin();
    else
      return simplexPool.begin(delta,stride);
  }
  
  const_simplexPtr_iterator simplexEnd(int delta=0, int stride=1) const
  {
    return simplexPool.end();
  }
  /*
  const_segmentPtr_iterator segmentBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return segmentPool.begin();
    else
      return segmentPool.begin(delta,stride);
  }
  
  const_segmentPtr_iterator segmentEnd(int delta=0, int stride=1) const
  {
    return segmentPool.end();
  }
  */
  const_vertexPtr_iterator vertexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return vertexPool.begin();
    else
      return vertexPool.begin(delta,stride);
  }
  
  const_vertexPtr_iterator vertexEnd(int delta=0, int stride=1) const
  {
    return vertexPool.end();
  }

  const_shadowSimplexPtr_iterator shadowSimplexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return shadowSimplexPool.begin();
    else
      return shadowSimplexPool.begin(delta,stride);
  }
  
  const_shadowSimplexPtr_iterator shadowSimplexEnd(int delta=0, int stride=1) const
  {
    return shadowSimplexPool.end();
  }

  const_shadowVertexPtr_iterator shadowVertexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return shadowVertexPool.begin();
    else
      return shadowVertexPool.begin(delta,stride);
  }
  
  const_shadowVertexPtr_iterator shadowVertexEnd(int delta=0, int stride=1) const
  {
    return shadowVertexPool.end();
  }

  const_ghostVertexPtr_iterator ghostVertexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return ghostVertexPool.begin();
    else
      return ghostVertexPool.begin(delta,stride);
  }
  
  const_ghostVertexPtr_iterator ghostVertexEnd(int delta=0, int stride=1) const
  {
    return ghostVertexPool.end();
  }

  const_ghostSimplexPtr_iterator ghostSimplexBegin(int delta=0, int stride=1) const
  {
    if (stride<2)
      return ghostSimplexPool.begin();
    else
      return ghostSimplexPool.begin(delta,stride);
  }
  
  const_ghostSimplexPtr_iterator ghostSimplexEnd(int delta=0, int stride=1) const
  {
    return ghostSimplexPool.end();
  }

  template <class L>
  void printIteratorWastedRatio()
  {
    if (glb::console->willPrint<L>())
      {
	glb::console->print<L>("Iterator wasted ratio :\n");
	glb::console->print<L>("  Vertices (l/g/s)  : %.2f / %.2f / %.2f\n",
			       vertexPool.getIteratorWastedRatio(),
			       ghostVertexPool.getIteratorWastedRatio(),
			       shadowVertexPool.getIteratorWastedRatio());
	glb::console->print<L>("  Simplices (l/g/s) : %.2f / %.2f / %.2f\n",
			       simplexPool.getIteratorWastedRatio(),
			       ghostSimplexPool.getIteratorWastedRatio(),
			       shadowSimplexPool.getIteratorWastedRatio());
      }
  }

protected:
  
  /** \brief Retrieve all simplices incident to the local vertices.
   * The simplices incident to vertex 'v' are stored in 'simplices', at indices
   * { index[v->getLocalIndex()] ... index[v->getLocalIndex()+1]-1}.
   * The number of simplices incident to vertex 'v' is therefore 
   * NS=index[v->getLocalIndex()+1]-index[v->getLocalIndex()]
   *  \param[out] simplices A vector that stores pointers to the simplices
   *  \param[out] index The index in 'simplices' of to the first simplex 
   *     incident to the nth vertex (i.e. vertex->getLocalIndex()==n)
   *  \param includeGhosts include ghost simplices if true, ignore them if false 
   *  \param nThreads number of threads to use.
   * \warning having more than one threads consumes much more memory !
   * \warning Be carefull that IT is large enough to store indices !
   */
  template <class IT>
  void getIncidentSimplices(std::vector<Simplex*> &simplices, 
			    std::vector<IT> &index,
			    bool includeGhosts=false,
			    int nThreads=glb::num_omp_threads)
  {
    // TimerPool::Timer timer;
    // timer.start();
    
    if (nThreads>1) 
      {
	getIncidentSimplices_parallel(simplices,index,includeGhosts,nThreads);
	return;
      }
    
    const unsigned long nVert=getNVertices();
    index.assign(nVert+1,0);
    IT * /*__restrict*/ indexP1 = &index[1];

    //glb::console->print<LOG_STD>("(%.3f)",timer.check());
    const simplexPtr_iterator it_end=simplexEnd();
    for (simplexPtr_iterator it=simplexBegin();it!=it_end;++it)
      {	
	IT id[Simplex::NVERT];
	
	for (int i=0;i<Simplex::NVERT;++i)
	  id[i]=it->getVertex(i)->getLocalIndex();

	for (int i=0;i<Simplex::NVERT;++i)
	  ++indexP1[id[i]];
      } 

    const ghostSimplexPtr_iterator git_end=ghostSimplexEnd();
    if (includeGhosts)
      {	

	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Not implemented correctly for ghost simplices.\n");
	glb::console->print<LOG_ERROR>("Need a fix!\n");
	exit(-1);
	//FIXME: THIS IS WRONG ! ghost vertices have different 
	for (ghostSimplexPtr_iterator git=ghostSimplexBegin();git!=git_end;++git)
	  {
	    for (int i=0;i<Simplex::NVERT;i++)
	      if (git->getVertex(i)->isLocal())
		{
		  index[git->getVertex(i)->getLocalIndex()+1]++;		    
		}	
	  } 
      }
    //glb::console->print<LOG_STD>("(%.3f)",timer.check());
    for (unsigned long i=1;i<=nVert;++i)
      index[i]+=index[i-1];
     
    simplices.resize(index.back());

    //glb::console->print<LOG_STD>("(%.3f)",timer.check());
   
    for (simplexPtr_iterator it=simplexBegin();it!=it_end;++it)
      {	
	
	IT id[Simplex::NVERT];
	for (int i=0;i<Simplex::NVERT;i++)
	  id[i]=it->getVertex(i)->getLocalIndex();
	for (int i=0;i<Simplex::NVERT;i++)
	  id[i]=index[id[i]]++;
	for (int i=0;i<Simplex::NVERT;i++)
	  simplices[id[i]]=(*it);
      } 

    if (includeGhosts)
      {
	for (ghostSimplexPtr_iterator git=ghostSimplexBegin();git!=git_end;++git)
	  {
	    for (int i=0;i<Simplex::NVERT;i++)
	      if (git->getVertex(i)->isLocal())
		{
		  simplices[index[git->getVertex(i)->getLocalIndex()]++]=(*git);
		}	
	  } 
      }

    //glb::console->print<LOG_STD>("(%.3f)",timer.check());
    for (unsigned long i=nVert;i>0;--i)
      index[i]=index[i-1];
    index[0]=0;
    //glb::console->print<LOG_STD>("(%.3f)",timer.check());
  }

  // remove stupid warning on GCC (gcc 4.4.7 prevent diagnostic inside functions ...
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

  template <class ICT>
  void updateIncidentSimplices(ICT &incidentSimplices,			       
			       bool includeGhosts=false,
			       int nThreads=glb::num_omp_threads)
  {
    typedef typename ICT::IT IT;
    std::vector<Simplex*> &simplices=incidentSimplices.simplices;
    std::vector<IT> &index=incidentSimplices.index;
    const unsigned long nVert=getNVertices();
    
    typedef unsigned short UpdtStatus;
    std::vector<UpdtStatus> updtStatus(nVert,0);    
    std::vector<Simplex*> newSimplices;
    std::vector<IT> newIndex(1,0);
    std::vector<IT> vertexIndex; // Stores the index of the updated vertex
  
#pragma omp parallel num_threads(nThreads)
    {
      int th = omp_get_thread_num();
      std::set<Simplex*> ball;
      std::vector<Simplex*> front;
      std::vector<Simplex*> newFront;
      
      // Variables for storing updated vertices incident simplices
      std::vector<Simplex*> newSimplicesLocal;
      std::vector<IT> newIndexLocal(1,0);
      std::vector<IT> vertexIndexLocal; // Stores the index of the updated vertex
      
      const simplexPtr_iterator it_end=simplexEnd();
      for (simplexPtr_iterator it=simplexBegin(th,nThreads);it!=it_end;++it)
	{
	  Simplex *s=(*it);
	  if (s->getLocalIndex() >= incidentSimplices.oldSimplicesCount)
	    {
	      for (int i=0;i<Simplex::NVERT;++i)
		{
		  Vertex *v=s->getVertex(i);
		  IT id=v->getLocalIndex();
		  // Check if this vertex is already taken care of
		  if (updtStatus[id]==0)
		    {
		      // NO -> check again with atomic operation for thread safety
		      unsigned char result;
		      //unsigned char *updtStatusPtr=&updtStatus[id];
		      //#pragma GCC diagnostic push
		      //#pragma GCC diagnostic ignored "-Wunused-value"

		      
		      /*
		      // Poor man's fetch and add ... for old intel compilers ...
#pragma omp critical
		      {			
			result=updtStatus[id];
			updtStatus[id] = updtStatus[id]+1;	
		      }
		      */
		      
#pragma omp atomic capture
		      { 
			// fetch_and_add
			result=updtStatus[id];
			updtStatus[id] = updtStatus[id]+1;			
		      }
		      //#pragma GCC diagnostic pop
		      
		      if (result==0)
			{
			  // NO conflict or we have priority
			  // identify the new ball around vertex id
			  ball.clear();
			  front.clear();

			  if (id<incidentSimplices.oldVerticesCount)
			    {	
			      for (long j=index[id];j<index[id+1];++j)
				if (simplices[j]->getVertexIndex(v)>=0)
				  front.push_back(simplices[j]);
			      
			      ball.insert(front.begin(),front.end());
			    }

			  if (ball.insert(s).second)
			    front.push_back(s);			
			  
			  do {
			    newFront.clear();
			    for (auto it=front.begin();it!=front.end();++it)
			      {
				for (int j=0;j<Simplex::NNEI;++j)
				  {
				    Simplex* nei=(*it)->getNeighbor(j);
				    
				    if ((nei!=NULL)&&
					(nei->isLocal()||includeGhosts)&&
					((*it)->getVertex(j)!=v)&&
					(ball.insert(nei).second))
				      newFront.push_back(nei);				
				  }
			      }
			    front.swap(newFront);
			  } while (front.size()>0);
			  
			  //printf("done.\n");
			  //and store the result for later use
			  vertexIndexLocal.push_back(id);			  
			  newIndexLocal.push_back(ball.size()+newIndexLocal.back());
			  newSimplicesLocal.insert(newSimplicesLocal.end(),
						   ball.begin(),ball.end());
			}
		    }
		}
	    }
	} // simplex iterator
      
      if (vertexIndexLocal.size()>0)
	{
#pragma omp critical
	  {	    
	    vertexIndex.insert(vertexIndex.end(),
			       vertexIndexLocal.begin(),vertexIndexLocal.end());
	    
	    newSimplices.insert(newSimplices.end(),
				newSimplicesLocal.begin(),newSimplicesLocal.end());
	    
	    long delta = newIndex.back();
	    long startIndex=newIndex.size();
	     newIndex.insert(newIndex.end(),
			     newIndexLocal.begin()+1,
			     newIndexLocal.end());
	   
	    for (long i=startIndex;i<newIndex.size();++i) 
	      newIndex[i]+=delta;	    
	  }
	}
    } // pragma omp
    
    index.resize(nVert+1);
    
    // update index for newly created vertices only.
    // After that, index will contain the index of the new adjacent simplices in 
    // newIndex and newSimplices
    long nNewVertexSimplices=0;
    long nUpdatedVertexSimplices=0;
  
    for (long i=0;i<vertexIndex.size();++i)
      {
	if (vertexIndex[i]>=incidentSimplices.oldVerticesCount)
	  {
	    // Newly created vertex
	    index[vertexIndex[i]+1]=i;
	    nNewVertexSimplices+=(newIndex[i+1]-newIndex[i]);
	  }
	else
	  {	   
	    // FIXME: There is a very rare bug were nNew is negative (=-1) causing an error
	    // have to fix this before reenabling upate in mesh.hxx ...

	    // updtStatus now contains the difference in the number of adjacent simplices
	    // before and after update
	    unsigned int nNew=
	      (newIndex[i+1]-newIndex[i])-
	      (index[vertexIndex[i]+1]-index[vertexIndex[i]]);
	    
	    static const long nNewMax=1L<<(sizeof(UpdtStatus)*8);
	    // It would be very suspicious if that happened ...
	    if (nNew>=nNewMax) {
	      PRINT_SRC_INFO(LOG_ERROR);
	      glb::console->print<LOG_ERROR>
		("More than %ld new incident simplices for one vertex !",nNewMax);
	      glb::console->print<LOG_ERROR>
		("This is highly suspicious, you may have a refinement problem ...");
	      glb::console->print<LOG_ERROR>
		("You can also change the type of 'updtStatus' to something wider if you think you understand what you're doing.\n");
	      exit(-1);
	    }
	    updtStatus[vertexIndex[i]]=1+nNew;
	    
	    nUpdatedVertexSimplices+=nNew;
	  }
      }
  
    long nNewSimplices=nNewVertexSimplices+nUpdatedVertexSimplices;
    simplices.resize(simplices.size()+nNewSimplices);
        
    //now we still have to take care of the updated vertices   
    //index[incidentSimplices.oldVerticesCount]-=nUpdatedVertexSimplices;
    long delta=nUpdatedVertexSimplices;
    std::vector<Simplex**> start;
    std::vector<Simplex**> stop;
    std::vector<Simplex**> out;
    start.reserve(nUpdatedVertexSimplices+1);
    stop.reserve(nUpdatedVertexSimplices+1);
    for (long i=incidentSimplices.oldVerticesCount-1;delta>0;--i)    
      {	
	if (updtStatus[i]==0)
	  {	  
	    // Store the range that needs to be copied
	    stop.push_back(&simplices[index[i+1]]);
	    start.push_back(&simplices[index[i+1]]);
	    out.push_back(stop.back()+delta);
	      
	    for (;updtStatus[i]==0;--i)
	      {
		start.back() -= (index[i+1]-index[i]);
		index[i+1]+=delta;	
	      }	  
	  }
	// This vertex was updated !
	// Leave free space to copy its simplices later and update index
	index[i+1]+=delta;
	delta-=(updtStatus[i]-1);
      }  

    int localNThreads = std::max(nThreads,2);
    // Remove spurious warning when OpenMP is not used
    UNUSED_VARIABLE(localNThreads);
#pragma omp parallel num_threads(localNThreads)
    {
#pragma omp single
      {
	// Copy the adjacent simplices for the new vertices
	// note that adjacent simplices are already put at their final location in simplices
#pragma omp task
	for (long i=incidentSimplices.oldVerticesCount;i<(index.size()-1);++i)
	  {
	    int id=index[i+1];
	    std::copy(&newSimplices[newIndex[id]],
		      &newSimplices[newIndex[id+1]],
		      &simplices[index[i]]);
	    index[i+1]=index[i]+newIndex[id+1]-newIndex[id];
	  }

#pragma omp task
	//#pragma omp parallel for num_threads(nThreads-1)
	  for (long i=0;i<start.size();++i)
	    std::copy_backward(start[i],stop[i],out[i]);
      }
    }
        
    // copy updated vertices adjacent simplices
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<vertexIndex.size();++i)
      {
	if (vertexIndex[i]<incidentSimplices.oldVerticesCount)
	  {
	    std::copy(&newSimplices[newIndex[i]],
		      &newSimplices[newIndex[i+1]],
		      &simplices[index[vertexIndex[i]]]);
	  }
      }         
  }

private:


  
  template <typename IT>
  void getIncidentSimplices_parallel(std::vector<Simplex*> &simplices, 
				     std::vector<IT> &index,
				     bool includeGhosts,
				     int nThreads)
  {
    TimerPool::Timer timer;
    timer.start();

    const unsigned long nVert=getNVertices();
    const unsigned long nIndices=nVert+1;
    std::vector< std::vector<IT> > thIndex(nThreads);

    index.resize(nIndices);
 
#pragma omp parallel num_threads(nThreads)
    {
      int th = omp_get_thread_num();
      
      thIndex[th].assign(nIndices,0);
      IT * /*__restrict*/ indexP1 = &thIndex[th][1];

      const simplexPtr_iterator it_end=simplexEnd();
      for (simplexPtr_iterator it=simplexBegin(th,nThreads);it!=it_end;++it)
	{
	  IT id[Simplex::NVERT];
	
	  for (int i=0;i<Simplex::NVERT;++i)
	    id[i]=it->getVertex(i)->getLocalIndex();

	  for (int i=0;i<Simplex::NVERT;++i)
	    ++indexP1[id[i]];
	}
      
      
      if (includeGhosts)
	{	
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Not implemented correctly for ghost simplices.\n");
	  glb::console->print<LOG_ERROR>("Need a fix!\n");
	  exit(-1);
	  //FIXME: THIS IS WRONG ! ghost vertices have different 
	  const ghostSimplexPtr_iterator git_end=ghostSimplexEnd();
	  for (ghostSimplexPtr_iterator git=ghostSimplexBegin(th,nThreads);
	       git!=git_end;++git)
	    {
	      for (int i=0;i<Simplex::NVERT;i++)
		if (git->getVertex(i)->isLocal())
		  {
		    indexP1[git->getVertex(i)->getLocalIndex()]++;	
		  }	
	    } 
	}

#pragma omp barrier
      unsigned long start=(nIndices/nThreads)*th;
      unsigned long stop=(nIndices/nThreads)*(th+1);
      if (th==nThreads-1) stop=nIndices;

      for (unsigned long i=start;i<stop;++i)
	index[i]=thIndex[0][i];
	
      for (unsigned long j=1;j<nThreads;++j)
	for (unsigned long i=start;i<stop;++i)
	  index[i]+=thIndex[j][i];
	
      for (unsigned long i=start;i<stop-1;++i)
	index[i+1]+=index[i];
      
#pragma omp barrier
      IT delta=0;
      for (unsigned long i=0;i<th;++i)
	delta+=index[(nIndices/nThreads)*(i+1)-1];

#pragma omp barrier
      for (unsigned long i=start;i<stop;++i)
	index[i]+=delta;

#pragma omp barrier
#pragma omp single
      {
	simplices.resize(index.back());
	//glb::console->print<LOG_STD>("(%.3f)",timer.check());
      }

      for (long j=nThreads-2;j>=0;--j)
	for (unsigned long i=start;i<stop;++i)
	  thIndex[j+1][i]=thIndex[j][i];
      //std::fill_n(&thIndex[0][start],stop-start,0);
      
      for (unsigned long j=1;j<nThreads-1;++j)
	for (unsigned long i=start;i<stop;++i)
	  thIndex[j+1][i] += thIndex[j][i];      
      
      if (start==0) start=1;
      std::copy_n(&index[start-1],stop-start,&thIndex[0][start]);

      for (unsigned long j=1;j<nThreads;++j)
	for (unsigned long i=start;i<stop;++i)
	    thIndex[j][i] += index[i-1];   
   
#pragma omp barrier
      //const simplexPtr_iterator it_end=simplexEnd();
      for (simplexPtr_iterator it=simplexBegin(th,nThreads);it!=it_end;++it)
	{
	  IT id[Simplex::NVERT];
	  for (int i=0;i<Simplex::NVERT;i++)
	    id[i]=it->getVertex(i)->getLocalIndex();
	  for (int i=0;i<Simplex::NVERT;i++)
	    {
	      id[i]=thIndex[th][id[i]+1]++;
	      //id[i]=indexP1[id[i]]++;
	    }
	  for (int i=0;i<Simplex::NVERT;i++)
	    simplices[id[i]]=(*it);
	  
	  //simplices[index[it->getVertex(i)->getLocalIndex()]++]=(*it);
	} 
      if (includeGhosts)
	{
	  const ghostSimplexPtr_iterator git_end=ghostSimplexEnd();
	  for (ghostSimplexPtr_iterator git=ghostSimplexBegin(th,nThreads);
	       git!=git_end;++git)
	    {
	      for (int i=0;i<Simplex::NVERT;i++)
		if (git->getVertex(i)->isLocal())
		  {
		    simplices[thIndex[th][git->getVertex(i)->getLocalIndex()+1]++]=(*git);
		  }	
	    } 
	}
    }
  
  }
  
};

/** \}*/
#include "../internal/namespace.footer"
#endif
