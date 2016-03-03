#ifndef __LOCAL_AMR_GRID_PROJECTOR_HXX__
#define __LOCAL_AMR_GRID_PROJECTOR_HXX__

#include <unistd.h>
#include <typeinfo>
#include <unordered_set>
#include <limits>

#include "../dice_globals.hxx"

#include "../tools/OMP/openMP_interface.hxx"

#include "../tools/helpers/helpers.hxx"

#include "../AMR/localAmrGridRaytracer.hxx"

#include "./internal/localAmrGridProjectorBase_2D.hxx"
#include "./internal/localAmrGridProjectorBase_3D.hxx"
#include "./internal/projectorWeightFunctor.hxx"
#include "./internal/localAmrGridContribSumInterface.hxx"
#include "./internal/concurrentQueue.hxx"
#include "../tools/helpers/fpuRoundingModeGuard.hxx"

/**
 * @file 
 * @brief  A class to project a simplicial mesh onto a local AMR grid. The mesh may be folded / auto-intersect, but
 *  should entirely fit within the grid boundaing box.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

//#define DEBUGDUMP 1
//7971373578124785664(0.499023,0.508789,0.864258)
/*
#define DEBUG_CHECK 7971373578124785664
#define DEBUG_CHECK_COORD {0.499023,0.508789,0.864258}
#define DEBUG_CHECK_RANK 0
*/

/*
//dbg2
#define DEBUG_CHECK 3701960856499648512
#define DEBUG_CHECK_COORD {0.499023,0.446289,0.401367}
//#define DEBUG_CHECK 5467372357105481728
//#define DEBUG_CHECK_COORD {0.499023,0.547852,0.592773}
#define DEBUG_CHECK_RANK 0
*/

//dbg3
/*
#define DEBUG_CHECK 9214367221307930624
#define DEBUG_CHECK_COORD {0.499023,0.541992,0.999023}
#define DEBUG_CHECK_RANK 0
*/
/*
#define DEBUG_CHECK 4629700425528045568
#define DEBUG_CHECK_COORD {0.591797,0.00195312,0.501953}
#define DEBUG_CHECK_RANK 0
*/


// 374
/*
#define DEBUG_CHECK 5665530354162726912
#define DEBUG_CHECK_COORD {0.499023,0.459961,0.614258}
#define DEBUG_CHECK_RANK 0
*/

// 403
/*
#define DEBUG_CHECK 5899717509016193024
#define DEBUG_CHECK_COORD {0.500977,0.454102,0.639648}
#define DEBUG_CHECK_RANK 0
*/
/*
#define DEBUG_CHECK 2283325012152680448
#define DEBUG_CHECK_COORD {0.500977,0.495117}
#define DEBUG_CHECK_RANK 0
#define DEBUG_CHECK_SINGLE_THREAD
*/
#ifdef DEBUGDUMP
#define PRINTVERTEX(v) {const Coord *crd=v->getCoordsConstPtr();printf("%e %e %e\n",crd[0],crd[1],crd[2]);}
#define PRINTCOORD(c) {printf("%e %e %e\n",c[0],c[1],c[2]);}
#endif

#ifdef DEBUG_CHECK
template <class VC>
class VoxelContribDebugger : public VC
{
public:
  void push_back(const typename VC::value_type& val)
  {
    
      static int count=0;
      long vid=DEBUG_CHECK;
      
      if ((val.first->index==vid)&&
	  (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))
	{
#pragma omp critical
	  {
	    std::cout.precision(20);
	    std::cout<<std::scientific<<++count<<": Adding "<<val.second
		     <<" to voxels "<<vid<<std::endl;
	  }
	}    
      VC::push_back(val);
  }
};
#else
template <class VC> class VoxelContribDebugger : public VC {};
#endif


/** \addtogroup AMR
 *   \{
 */

/**
 * \class ProjectionTag
 * \brief A policy class used to set the projection mode of the simplices
 */
class ProjectionTag {
public:
  enum Type {
    regular=0,/*!< Use exact integration */
    sampleNoOverlap=(1<<0), /*!< Project the simplex as though it was a single point */
    sampleWithOverlap=(1<<1), /*!< Project the simplex using sampling points */
    sampleFromVertices=(1<<2), /*!< Use vertices as mass tracers to project the simplex */
    skip=(1<<3) /*!< ignore the simplex for projection */
  };
  
  ProjectionTag(double maxLevel=-1, Type lvlMethod=sampleFromVertices, 
		double minVolume=-1, Type vMethod=sampleWithOverlap, 
		double minLength=-1, Type lMethod=sampleWithOverlap, 
		double maxAnisotropy=-1, Type aMethod=sampleWithOverlap)
  {
    this->levelThreshold=
      (maxLevel<0)?(std::numeric_limits<double>::max()):maxLevel;
    this->anisotropyThreshold=
      (maxAnisotropy<0)?(std::numeric_limits<double>::max()):maxAnisotropy;
    this->volumeThreshold=minVolume;
    this->lengthThreshold=minLength;

    this->lvlMethod=lvlMethod;
    this->aMethod=aMethod;
    this->vMethod=vMethod;
    this->lMethod=lMethod;
  }
  
  template <class AMR, class Simplex>
  void tagSimplex(AMR *amr,Simplex *s, 
		  double bBox[2][Simplex::NDIM],
		  double maxLevelVoxelLengthInv[Simplex::NDIM],
		  double v, double sMin, double sMax) const
  {
    s->cache.c[0]=ProjectionTag::regular;
    bool inside=true;
    for (int j=0;j<Simplex::NDIM;++j)
      {		  
	int a=(int)((bBox[0][j]-amr->getBBoxMin(j)-amr->getEpsilon(j))*
		    maxLevelVoxelLengthInv[j]);
	int b=(int)((bBox[1][j]-amr->getBBoxMin(j)+amr->getEpsilon(j))*
		    maxLevelVoxelLengthInv[j]);
	inside &= (a==b);
      }
   
    if (inside) s->cache.c[0]=ProjectionTag::sampleNoOverlap;    
    else 
      {
	if (s->getLevel()>=levelThreshold) //level
	  s->cache.c[0]=lvlMethod;	
	else if (v<volumeThreshold) //volume
	  s->cache.c[0]=vMethod;
	else if (v<sMax*lengthThreshold)  // minimum width
	  s->cache.c[0]=lMethod;
	else if (sMax > sMin*anisotropyThreshold) //anisotropy
	  s->cache.c[0]=aMethod;
      }    
  }
  /*
  template <class Simplex>
  static void tagSimplex(Simplex *s, Type tag)
  {
    s->cache.c[0]=tag;
  }
  */
private:
  double anisotropyThreshold;
  double volumeThreshold;
  double lengthThreshold;
  double levelThreshold;

  Type aMethod;
  Type vMethod;
  Type lMethod;
  Type lvlMethod;
};

/**
 * \class LocalAmrGridProjectorT
 * \brief A class to project a mesh onto a local AMR grid. The unstructured mesh MUST be
 * contained within the bounds of the AMR grid (always true for periodic boundaries)
 * We use simulation of simplicity to deal with degeneracies, with the convention that the
 * voxel's vertices ith coordinate is perturbed as follows:
 *     Xi -> Xi - pow(eps,NDIM+1-i)
 * \tparam AMR  The local AMR grid class
 * \tparam MESH The mesh class
 * \tparam checkAccuracy Enable/Disable accuracy checking. Note that checking accuracy may
 * limit the number of threads you can use efficiently in parallel (a warning will be 
 * issued if this happend).
 * \tparam IF   A floating point type to use for the computation of the individual 
 * contributions to each voxel. Defaults to 'long double'.
 * \tparam HF   A high precision floating point type to use for the computation of 
 * the individual contributions to each voxel when IF is not precise enough. Defaults 
 * to 'long double'.
 * \tparam SF   A floating point type to use for the summation of all the individual 
 * contributions to each voxel. If this is different from Amr::Data, then an additional 
 * temporary array of size sizeof(SF)*number_of_voxels will have to be allocated. Defaults 
 * to 'long double'.
 */

template <class AMR, class MESH, bool checkAccuracy=false,
	  class IF=long double, 
	  class HF=long double, 
	  class SF=long double> 
class LocalAmrGridProjectorT : 
  protected internal::LocalAmrGridProjectorBaseT<AMR::NDIM,AMR,MESH,IF,HF,checkAccuracy,0>,
  protected internal::LocalAmrGridProjectorBaseT<AMR::NDIM,AMR,MESH,HF,HF,checkAccuracy,1>
{
private:
  typedef internal::LocalAmrGridProjectorBaseT<AMR::NDIM,AMR,MESH,IF,HF,checkAccuracy,0> 
  IFBase;
  typedef internal::LocalAmrGridProjectorBaseT<AMR::NDIM,AMR,MESH,HF,HF,checkAccuracy,1> 
  HFBase;

  typedef internal::LocalAmrGridContribSumInterfaceT<AMR,SF,HF,checkAccuracy> 
  ContribSumInterface;

  //typedef typename Base::IncidentEdge IncidentEdge;  

public:    
  static const int NDIM = AMR::NDIM; 

  //typedef IF IFloat;
  //typedef HF HFloat;

  typedef typename MESH::vertexPtr_iterator vertexPtr_iterator;
  typedef typename MESH::simplexPtr_iterator simplexPtr_iterator;
  typedef typename MESH::Simplex Simplex;
  typedef typename MESH::Vertex  Vertex;
  typedef typename MESH::Segment Segment;
  typedef typename MESH::SegmentHandle SegmentHandle;
  typedef typename MESH::Facet Facet;
  typedef typename MESH::FacetHandle FacetHandle;
  typedef typename MESH::Coord Coord;
  typedef typename MESH::GeometricProperties MeshGeometricProperties;    

  typedef typename AMR::ICoord ICoord;
  typedef typename AMR::Voxel Voxel;
  typedef typename AMR::GeometricProperties AmrGeometricProperties;  
  typedef typename AMR::Data AmrData;

  typedef MeshGeometricProperties MGP;
  typedef AmrGeometricProperties AGP;

  //typedef LocalAmrGridRaytracerT<AMR,IF> Raytracer;
  //typedef typename Raytracer::Coord RayCoord;

  /**
   *  \param amr_ A pointer to the AMR grid 
   *  \param mesh_ A pointer to the unstructured mesh to project
   *  \param nThreads_ the number of threads to use
   *  \param verbose_ whether to output progress on the console
   *  \pre The mesh must be contained within the bounds of the amr grid.
   */
  LocalAmrGridProjectorT(AMR *amr_, MESH *mesh_, 
			 double accLevel=checkAccuracy?0.1:0,
			 int nThreads_=glb::num_omp_threads,
			 bool verbose_=false):			 
    IFBase(amr_,mesh_,nThreads_),
    HFBase(amr_,mesh_,nThreads_),
    amr(amr_),
    mesh(mesh_),
    nThreads(nThreads_),
    //mpiCom(mpiCom_),
    verbose(verbose_),
    samplesPerVoxel(10),
    contribSumInterface(amr,accLevel)
    //highPrecisionBase(amr_,mesh_,nThreads_)
  {
#ifdef DEBUG_CHECK_SINGLE_THREAD
    nThreads=1;
#endif


    amrGeometry = amr->getGeometry();
    meshGeometry = mesh->getGeometry();
    
    int tmpSegDir[Voxel::NSEGNEI][NDIM];
    // Precompute the neighborDirection table
    for (int i=0;i<AMR::Voxel::NVERT;++i)
      {
	amr->getVertexNeighborsDirection(i,vertexNeighborDir[i]);
	for (int j=0;j<NDIM;++j)
	  {
	    amr->getSegmentNeighborsDirection(i,j,tmpSegDir);
	    for (int k=0;k<NDIM;++k)
	      segmentNeighborDir[i][j][k] = tmpSegDir[0][k];
	  }
      }
    
    if (NDIM==2) 
      dimFactor  = 0.5;
    else if (NDIM==3) 
      dimFactor = -1.0/6.0;

    initProjectorTimer = glb::timerPool->pop("amrProjector_initProjection");
    initTimer = glb::timerPool->pop("amrProjector_init");
    projectionTimer = glb::timerPool->pop("amrProjector_project");
    reprojectionTimer = glb::timerPool->pop("amrProjector_reproject");

#ifdef DEBUGDUMP
    flct=0;
    nDuplicates=0;
    openFile("project");
    openFile("edges");
    openFile("voxelVertices",true);
#endif
    /*
    static int ct=0;
    char name[255];
    sprintf(name,"particles_%5.5d.dat",ct++);
    pf=fopen(name,"w");
    */

    //countT=0;
    /*
    for (int i=0;i<AMR::Voxel::NVERT;++i)
      {
	printf("VERTEX %d:\n",i);
	for (int j=0;j<(1<<NDIM);++j)
	  printf("(%d,%d)",vertexNeighborDir[i][j][0],vertexNeighborDir[i][j][1]);
	printf("\n");	    
      }
    exit(0);
    */
  }
  //FILE *pf;
  ~LocalAmrGridProjectorT()
  {
    //fclose(pf);
#ifdef DEBUGDUMP
    closeFiles();
#endif
  }

  /** \brief Advancement of projection will be printed to conole if state is set to true
   */
  void setVerbose(bool state)
  {
    verbose=state;
  }
 
  /** 
   * \brief project the mesh onto the AMR grid with the weight function approximated
   * to order 0 (i.e. constant over each simplex)  
   * \param wf a functor that returns the average weight over a SIMPLEX (see WF)
   * \param checkTags if true, read the simplices cache.c[0] value and interpret it as 
   * follows: 0: nothing special, 1: simplex is fully contained inside a voxel
   * 2: simplex is close to degenerate
   * \return the number of simplices that had to be reprojected
   * \tparam WF the class of the weight functor that implements a 
   * WF::operator()(MESH::Simplex *s) method returning the average weight of a simplex \a s.
   * \remark The simplices cache is invalidated
   */
  template <class WF>
  long project(const WF &wf, bool checkTags=true)
  {
    long nReproj=0;
    glb::console->printFlush<LOG_INFO>("Initializing AMR grid projector ... ");
    initProjectorTimer->start();
    //tags are needed for checking accuracy !
    if ((!checkTags)&&(!checkAccuracy))
      {
	typedef internal::ProjectorWeightFunctorT
	  <0,MESH,WF,WF,IF,HF> WeightFunctor;
	WeightFunctor weightFunctor(mesh,wf,nThreads);
	glb::console->printFlush<LOG_INFO>("done in %lgs.\n",initProjectorTimer->stop());
	nReproj=run<WeightFunctor>(weightFunctor);	
      }
    else
      {
	// If we correct for lack of accuracy, we'll need the full accuracy for the
	// weight also ...
	typedef typename hlp::IF_<checkAccuracy,HF,IF>::Result IFloat;

	typedef internal::ProjectorWeightFunctorT
	  <0,MESH,WF,WF,IFloat,HF,internal::Tag_ExcludeIfTagged> 
	  WeightFunctor;
	WeightFunctor weightFunctor(mesh,wf,nThreads,!checkTags);
	glb::console->printFlush<LOG_INFO>("done in %lgs.\n",initProjectorTimer->stop());
	nReproj=run<WeightFunctor>(weightFunctor);
      }
    return nReproj;
  }

  /** 
   * \brief project the mesh onto the AMR grid with the weight function approximated
   * to order 1 (i.e. linearly interpolated over each simplex)  
   * \param wf a functor that returns the value of the weight at VERTICES locations 
   *  (vertex is given as argument to the functor, which returns its weight, see WF)
   * \param wdf a functor that sets the derivative of the weight over a SIMPLEX (see WDF).
   * \return the number of simplices that had to be reprojected
   * \tparam WF The class of the weight functor that implements a 
   * WF::operator()(MESH::Vertex *v) method returning the weight evaluated at the location
   * of vertex \a v.
   * \tparam WDF the class of the weight derivative functor that implements a 
   * WDF::operator()(MESH::Simplex *s, double result[NDIM]) method that sets the NDIM values
   * of result to the derivatives of the weight within the simplex \a s.
   * \param checkTags if true, read the simplices cache.c[0] value and interpret it as 
   * follows: 0: nothing special, 1: simplex is fully contained inside a voxel
   * 2 or 3: simplex is very close to degenerate
   * \remark The simplices cache is invalidated
   */
  template <class WF, class WDF>
  long project(const WF &wf, const WDF &wdf, bool checkTags=true)
  {    
    long nReproj=0;   
    glb::console->printFlush<LOG_INFO>("Initializing AMR grid projector ... ");
    initProjectorTimer->start();
    //tags are needed for checking accuracy !
    if ((!checkTags)&&(!checkAccuracy))  
      {
	internal::ProjectorWeightFunctorT
	  <1,MESH,WF,WDF,IF,HF>
	  weightFunctor(mesh,wf,wdf,nThreads);
	glb::console->printFlush<LOG_INFO>("done in %lgs.\n",initProjectorTimer->stop());
	nReproj=run(weightFunctor);
      }
    else
      {
	// If we correct for lack of accuracy, we'll need the full accuracy for the
	// weight also ...
	typedef typename hlp::IF_<checkAccuracy,HF,IF>::Result IFloat;

	internal::ProjectorWeightFunctorT
	  <1,MESH,WF,WDF,IFloat,HF,internal::Tag_ExcludeIfTagged> 
	  weightFunctor(mesh,wf,wdf,nThreads);
	glb::console->printFlush<LOG_INFO>("done in %lgs.\n",initProjectorTimer->stop());
	nReproj=run(weightFunctor);
      }
    return nReproj;
  }

  /** 
   * \brief sets the number of sampling points used when projecting simplices with mode
   * ProjectionTag::sampleWithOverlap
   */
  void setSamplesPerVoxel(long nSamples)
  {
    samplesPerVoxel=nSamples;
  }

  /** 
   * \brief get the number of sampling points used when projecting simplices with mode
   * ProjectionTag::sampleWithOverlap
   */
  long getSamplesPerVoxel() const
  {
    return samplesPerVoxel;
  }

private:
  AMR *amr;
  MESH *mesh;
  int nThreads;
  //MpiCommunication *mpiCom;

  MeshGeometricProperties *meshGeometry;
  AmrGeometricProperties *amrGeometry;
  int vertexNeighborDir[AMR::Voxel::NVERT][AMR::Voxel::NVERTNEI][NDIM];
  int segmentNeighborDir[AMR::Voxel::NVERT][NDIM][NDIM];
  AmrData dimFactor;
  //double dimFactor1;
  //double dimFactorInv;
  //double dimFactorInv1;
  
  typename TimerPool::Timer *initProjectorTimer;
  typename TimerPool::Timer *initTimer;
  typename TimerPool::Timer *projectionTimer;
  typename TimerPool::Timer *reprojectionTimer;

  int busyFlag;
  internal::ConccurentQueue cQueue;
  bool verbose;
   
  long samplesPerVoxel;

#ifdef DEBUGDUMP
  // Commented because it generates warnings on GCC even when DEBUGDUMP is not defined
  // I don't understand why ...
  /*
  long nDuplicates;
  std::set< std::pair<typename Voxel::ICoord,int> > vRefSet;
  long foundVertices;
  int flct;
 
  FILE *fl[10];
  char fname[10][255];
  void openFile(const char *name, bool increaseCount=false)
  {
    static int count =0;
    //char fname[255];
    sprintf(fname[flct],"%s_%6.6d.dat",name,count);
    fl[flct]=fopen(fname[flct],"w");
    fprintf(fl[flct],"ANDNET\n3\n?????????\n");
    flct++;
    if (increaseCount) count++;
  }

  void closeFiles()
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtrigraphs"

    for (int i=0;i<flct;++i)
      {
	fclose(fl[i]);
	char cmd[256];
	sprintf(cmd,"cat %s | sed \"s/?????????/$(( `cat %s |wc -l` - 3))/g\" > %s.OK ",
		fname[i],fname[i],fname[i]);
	system (cmd);
	sprintf(cmd,"mv %s.OK %s",fname[i],fname[i]);
	system (cmd);
      }
    flct=0;
#pragma GCC diagnostic pop
  }  
  */
#endif

  /** 
   * \brief handle the initialization, threads, accuracy ... for the projection
   * \return the number of simplices that had to be reprojected
   * The actual projection is done in do_project.
   * \note The correct FPU rounding mode is automatically set if the floating point types
   * require it.
   */
  template <class WF>
  long run(const WF &weightFunctor)
  {   
    long nRepSimplices=0;
    int pass=0;
    int bAcc;
    double dAcc;
    bool logToStdOnly = ((!glb::console->willPrint<LOG_INFO>())&&
			 (glb::console->willPrint<LOG_STD>()));

    // set fpu rounding mode if necessary, depending on the types ...
    hlp::FpuRoundingModeGuardT<IF> fpuGuardgI;
    hlp::FpuRoundingModeGuardT<HF> fpuGuardgH;
    hlp::FpuRoundingModeGuardT<SF> fpuGuardgS;

    contribSumInterface.getAccuracyLevel(bAcc,dAcc);
    if (verbose)
      {
	if (logToStdOnly)
	  {
	    glb::console->printFlush<LOG_STD>
	      ("Projecting mesh at order %d, acc=%g(%d bits), s=%d/%d bits ... ",
	       WF::ORDER,dAcc,bAcc,
	       std::numeric_limits<IF>::digits,
	       std::numeric_limits<HF>::digits);
	  }
	else
	  {
	    glb::console->print<LOG_INFO>
	      ("Projecting mesh at order %d, acc=%g(%d bits), s=%d/%d bits:\n",
	       WF::ORDER,dAcc,bAcc,
	       std::numeric_limits<IF>::digits,
	       std::numeric_limits<HF>::digits);
	    //("Projecting mesh (acc=%g / %d bits, order %d):\n",dAcc,bAcc,WF::ORDER);
	    glb::console->indent();
	    glb::console->printFlush<LOG_INFO>("* Initializing mesh and grid ... ");
	  }
      }

    initTimer->start();
    //openFile("project");    
    
    // Retrieve the simplices incident to each vertex
    if (verbose) glb::console->printFlush<LOG_PEDANTIC>("(incidence)");     
   
    typedef typename Vertex::LocalIndex LocalIndex;

    auto incidence=mesh->getIncidentSimplices();
    const std::vector<Simplex*> &incidentSimplices=incidence.first;
    const std::vector<LocalIndex> &iIndex=incidence.second;
    
    double incidenceTime = initTimer->check();   

    // Identify individual vertices and segments of the AMR grid
    if (verbose) glb::console->printFlush<LOG_PEDANTIC>("(AMR)");
    // note : we only need voxel segments for NDIM==3 (hence the NDIM>2)
    amr->assignVerticesToLeaves(NDIM>2,nThreads);

    double initTime=initTimer->stop();
    if (verbose) 
      {
	glb::console->print<LOG_INFO>
	  (" done in %lgs. (%.1f/%.1f%%)\n",initTime,
	   100.0*incidenceTime/initTime,
	   100.0*(1.0-incidenceTime/initTime));
      }

    double projectionTime=0;
  reproject: // Reprojection starting point !

    projectionTimer->start();
#ifdef DEBUGDUMP
    foundVertices=0;
    glb::console->print<LOG_DEBUG>
      (" Counted %ld unique vertices, %ld unique segments.\n",
       amr->getUniqueVerticesCount(),
       amr->getUniqueSegmentsCount());
#endif

    double reprojectionInitTime=0;
    if (pass>0)
      {
	typename TimerPool::Timer timer;
	timer.start();

	glb::console->printFlush<LOG_INFO>("* Initializing reprojection ... ");
	
	// We are reprojecting failed voxels
	FOREACH_THREAD_SIMPLEX(mesh,nThreads,th,it)
	  {
	    Coord simplexBBox[2][NDIM];
	    std::vector<Voxel*> overlap;
	    long nRepSimplicesLocal=0;
	    
	    for (;it!=it_end;++it)
	      {
		Simplex *simplex=(*it);
		const Coord *refCoords=simplex->getVertex(0)->getCoordsConstPtr(); 
		for (int j=0;j<NDIM;++j) 
		  simplexBBox[0][j]=simplexBBox[1][j]=refCoords[j];
		for (int j=1;j<Simplex::NVERT;++j)
		  {
		    Vertex *v=simplex->getVertex(j);
		    const Coord *coords=v->getCoordsConstPtr();
		    for (int k=0;k<NDIM;++k)
		      {
			// We have to be a bit carefull with periodic boundaries !
			Coord c = meshGeometry->checkCoordConsistency
			  (coords[k],refCoords[k],k);
	    
			if (simplexBBox[0][k]>c)
			  simplexBBox[0][k]=c;
			if (simplexBBox[1][k]<c)
			  simplexBBox[1][k]=c;
		      }	
		  }

		overlap.clear();
		amr->getLeavesBBoxOverlap(simplexBBox,std::back_inserter(overlap));
		bool overlaps=false;
		for (unsigned long ovi=0;ovi<overlap.size();ovi++)
		  {
		    if (contribSumInterface.needReprojection(overlap[ovi]))
			overlaps=true;
		  }
		if (!overlaps) 
		  weightFunctor.setTag(simplex,ProjectionTag::skip);
		else
		  nRepSimplicesLocal++;
	      }
#pragma omp atomic
	    nRepSimplices+=nRepSimplicesLocal;
	  }
	FOREACH_THREAD_END;

	reprojectionInitTime = timer.stop();
	glb::console->printFlush<LOG_INFO>("done in %lgs. (%ld simplices tagged)\n",
					   reprojectionInitTime,nRepSimplices);
      }

    // Now do the actual projection ...
    long totalContribs = 0;
    if (pass==0)
      totalContribs=projectionHandler<IFBase,0>(incidentSimplices,iIndex,weightFunctor);
    else
      totalContribs=projectionHandler<HFBase,1>(incidentSimplices,iIndex,weightFunctor);

    // And check accuracy
    int nFailed = 0;
    if (pass==0)
      {
	reprojectionTimer->start();
	if (verbose) 
	  {
	    if (!logToStdOnly)	  
	      {	    
		double elapsed=projectionTimer->check();
		glb::console->print<LOG_INFO>("done in %lgs.\n",elapsed);	    
		glb::console->print<LOG_PEDANTIC>
		  ("* Total number of intersections: %ld.\n",totalContribs);
		glb::console->print<LOG_PEDANTIC>
		  ("* Intersections processed per thread per second : %lg/th/s.\n",
		   double(totalContribs)/elapsed/nThreads);	    
		glb::console->printFlush<LOG_INFO>("* Checking accuracy ... ");
	      }
	  }   

	nFailed = contribSumInterface.setReprojectionModeIfNeeded(nThreads);

	if (verbose) 
	  {
	    if (logToStdOnly)
	      {
		if (nFailed>0)
		  glb::console->print<LOG_STD>("+%d voxels ... ",nFailed);	    	 
	      }
	    else
	      {
		if (checkAccuracy)
		  glb::console->printFlush<LOG_INFO>("%d inaccurate voxels found in %gs.\n",
						     nFailed,reprojectionTimer->check());
		else
		  glb::console->printFlush<LOG_INFO>("SKIPPED.\n");
	      }
	  }
	// Reproject if we failed somewhere ...
	if (nFailed>0)
	  {	
	    if (hlp::SameType<IF,HF>::value)
	      {		
		glb::console->printFlush<LOG_STD>("\n");		
		glb::console->printFlush<LOG_WARNING>
		  ("Accuracy criteria for projection could not be reached and higher numerical precision is not available.\n");
	      }
	    else
	      {
		pass++;
		projectionTime=projectionTimer->stop();		
		goto reproject;
	      }
	  }
	reprojectionTimer->stop();
      }
    else
      {
	if (pass>0) reprojectionTimer->stop();
	if (verbose) 
	  {
	    if (!logToStdOnly)
	      {	    
		double elapsed=projectionTimer->check()-reprojectionInitTime;
		glb::console->print<LOG_INFO>("done in %lgs.\n",elapsed);	    
		glb::console->print<LOG_PEDANTIC>
		  ("* Total number of intersections: %ld.\n",totalContribs);
		glb::console->print<LOG_PEDANTIC>
		  ("* Intersections processed per thread per second : %lg/th/s.\n",
		   double(totalContribs)/(elapsed)/nThreads);
	      }
	  }   
      }

    // Assign the computed values to the voxels and apply correction factor
    contribSumInterface.commit(nThreads);
    amr->visitTree(CorrectDimFactorVisitorT<true>(amr,dimFactor),nThreads); 

    projectionTime+=projectionTimer->stop();  
    if (verbose) 
      {
	if (logToStdOnly)
	  {
	    glb::console->print<LOG_STD>("done. (%lgs / %lgs)\n",
					 initTime,projectionTime);	    	    
	  }
	else
	  {
	    glb::console->unIndent();
	    glb::console->print<LOG_INFO>("Projection done in %lgs.\n",
					  initTime+projectionTime);
	  }
      }    

#ifdef DEBUGDUMP    
    long N=64+1; // that's for a 64*64*64 only
    long NN = (NDIM==3)?(6*N*N-12*N+8):4*N;
    glb::console->print<LOG_STD_ALL>(" Found %ld(+%ld=%ld) vertices (%ld total)\n",
				     foundVertices,NN,NN+foundVertices,
				     amr->getUniqueVerticesCount());
#endif
    //closeFile();

    return nRepSimplices;
  }
  
  template <class Base, int Pass, class WF, typename IT>
  long projectionHandler(const std::vector<Simplex*> &incidentSimplices,
			 const std::vector<IT> &iIndex,			 
			 const WF &weightFunctor)
  {
    long contribCount[nThreads];
    if (nThreads>1)
      {
	// Max number of contribs to store per thread before commit
	// We allow 4Mb per thread
	const long nReserved = (1<<22)/sizeof(typename Base::Float);
	// We need at least 128 simplices per batch so that
	// multithreading is meaningfull ...
	long nBatches = mesh->getNVertices()/128; 
	if (nBatches < 16*nThreads) nBatches=16*nThreads;
	if (nBatches > 512*nThreads) nBatches=512*nThreads;

	if (verbose) 
	  glb::console->printFlush<LOG_INFO>
	    ("* Projecting (%d threads/%ld batches) ... ",nThreads,nBatches);

	// Debugger is empty if DEBUG_CHECK is not defined
	typedef VoxelContribDebugger<
	  std::vector< std::pair<Voxel*,typename Base::Float> > 
	  > VoxelContribs;    
	
	static const long voxelContribsStride=16;
	std::vector< VoxelContribs > voxelContribs_threads(nThreads*voxelContribsStride);

	for (int i=0;i<nThreads;++i) 
	  {
	    voxelContribs_threads[i*voxelContribsStride].reserve(nReserved);
	    contribCount[i]=0;
	  }
	busyFlag=0;

	cQueue.clear();
	cQueue.resize(nThreads);
	bool warn=false;
	
#pragma omp parallel num_threads(nThreads)
	{
	  int th=omp_get_thread_num();
	  VoxelContribs &voxelContribs=voxelContribs_threads[th*voxelContribsStride];
	  
#pragma omp for schedule(dynamic) nowait
	  for (long bId=0;bId<nBatches;++bId)
	    {
	      std::vector<Voxel*> overlap;
	   
	      vertexPtr_iterator it=mesh->vertexBegin(bId,nBatches);	    
	      const vertexPtr_iterator it_end=mesh->vertexEnd(bId,nBatches); 
	      do
		{	
		  // The maximum accumulated contributions to allow
		  long nMax = nReserved*1.9 - voxelContribs.size();
		  // The trigger to ask for flushing contributions
		  long enrolTrigger=0;
		  if (voxelContribs.size() < nReserved/2) 
		    enrolTrigger=(nReserved/2) - voxelContribs.size();
		
		  // Prevent the contribs from growing too much when contribSumInterface
		  // cannot keep up the pace because we are using too many threads
		  if (nMax>0)
		    do_project<Base,Pass>(incidentSimplices,iIndex,overlap,
					  weightFunctor,std::back_inserter(voxelContribs),
					  it,it_end,enrolTrigger,nMax,th);
		  else if (!warn)
		    {
		      warn=true;
		      glb::console->print<LOG_WARNING>("\n");
		      glb::console->print<LOG_WARNING>
			("Thread %d waiting: The (single threaded) summation cannot always keep up the pace with contribs computation !",th);
		      glb::console->print<LOG_WARNING>
			("You may be using too many threads for projection (%d). Consider lowering threads count if you observe important performance degradation ...\n",
			 nThreads);
		    }
		
		  if (cQueue.availableFor()==th)
		    {
		      //printf("Processing %d, %ld waiting.\n",th,cQueue.nWaiting());
		      const unsigned long contribsSize=voxelContribs.size();
		      for (unsigned long j=0;j<contribsSize;++j)
			contribSumInterface.add(voxelContribs[j]);

		      contribCount[th]+=contribsSize;
		      voxelContribs.clear();
		      cQueue.done(th);
		    }
		
		} while (it!=it_end);
	    }

	  cQueue.enrol(th);
          while (cQueue.isEnroled(th))
            {
              if (cQueue.availableFor()==th)
                {
                  //printf("Finishing %d, %ld waiting.\n",th,cQueue.nWaiting());
                  const unsigned long contribsSize=voxelContribs.size();
                  for (unsigned long j=0;j<contribsSize;++j)
                    contribSumInterface.add(voxelContribs[j]);

                  contribCount[th]+=contribsSize;
                  voxelContribs.clear();
                  cQueue.done(th);
                  //printf("Processed %d, next in line: %d\n",th,cQueue.availableFor());
                }
            }	 
	} // parallel region
      }
    else
      {
	// Single threaded version
	if (verbose) glb::console->printFlush<LOG_INFO>("* Projecting (%d thread) ... ",nThreads);
	vertexPtr_iterator it=mesh->vertexBegin();
	const vertexPtr_iterator it_end=mesh->vertexEnd();
	std::vector<Voxel*> overlap;

	FakeBackInserterT<ContribSumInterface> out(contribSumInterface);	
	contribCount[0]=do_project<Base,Pass>
	  (incidentSimplices,iIndex,overlap,weightFunctor,out,it,it_end);
      }   
    long totalContribs=0;
    for (int i=0;i<nThreads;++i) 
      totalContribs+=contribCount[i];

    return totalContribs;
  }

  // NOTE: This function is thread safe
  // NOTE: simplex vertices that fall on a right boundary of the AMR box should NOT 
  // contribute when boundary conditions are not periodic. However this is difficult to 
  // handle so we let them contribute, but compensate in addVoxelsContrib by removing
  // the contribution of the voxel's vertices falling on the right boundary.
  template <class Base, int Pass, class WF, class OutputIterator, typename IT>
  long do_project(const std::vector<Simplex*> &incidentSimplices,
		  const std::vector<IT> &iIndex,
		  std::vector<Voxel*> &overlap,
		  const WF &weightFunctor,
		  OutputIterator out, 
		  vertexPtr_iterator &it,
		  const vertexPtr_iterator &it_end,
		  long enrolTrigger=-1,
		  long nMax=-1,
		  int threadIndex=0)
  {
    typedef typename Base::Float  IFloat;
    typedef typename Base::HFloat HFloat;

    typedef LocalAmrGridRaytracerT<AMR,IFloat> Raytracer;
    typedef typename Raytracer::Coord RayCoord;
    typedef typename Base::IncidentEdge IncidentEdge;

    Raytracer raytracer(amr); // Define it here as project may be called within threads
    std::vector<IncidentEdge> ownedEdges;
    long nProcessed=0;
  
    /* // for debugging 
    long guessHPSuccess=0;
    long guessHPFailure=0;
    long guessHPTotalFailure=0;
    long nHPNeeded=0;
    long nHPComputed=0;
    std::cout.precision(15);
    std::cout << std::scientific;
    */	

#ifdef DEBUG_CHECK
#ifdef DEBUG_CHECK_COORD
#pragma omp critical
    if (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK)
      {	
	Coord c[NDIM]=DEBUG_CHECK_COORD;
	Voxel *v = amr->getVoxelAt(c);
	
	if (v->getIndex()!=DEBUG_CHECK) 
	  {
	    glb::console->print<LOG_STD_ALL>("Found voxel ID : %ld\n",v->getIndex());
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_STD_ALL>
	      ("Debug mode activated, but DEBUG_CHECK and DEBUG_CHECK_COORD are not consistent !\n");
	    glb::console->print<LOG_STD_ALL>
	      ("Update DEBUG_CHECK_COORD to the coordinates of a point inside the voxel with id DEBUG_CHECK or deactivate debug by undefining DEBUG_CHECK.\n");
	    exit(-1);
	  }	
      }
#endif
#endif
    
    
    // loop over all the vertices in the tesselation
    for (;it!=it_end;++it)
      {	
	// for debugging
	//int guessHP=false;
	//int needHP=false;
	
	if ((enrolTrigger>=0)&&(nProcessed>=enrolTrigger))
	  {	   
	    cQueue.enrol(threadIndex);
	    if ((cQueue.availableFor()==threadIndex)||(nProcessed>nMax))
	      break;	   
	  }
	
	/*
	// Check the number of contributions so far and break the loop if there are
	// too many, so that we can process them before we continue 
	if (nProcessed>maxCount) 
	  {
	    break;
	    // Add a random number of vertices to each process so that threads do not get
	    // accidentally synchronized ?
	    
	    // if (busyFlag) maxCount+=10+((maxCount+65537)%100);
	    // else break;
	    
	  }
	*/
	Vertex *vertex = (*it);
	Voxel *voxel = amr->getVoxelAt(vertex->getCoordsConstPtr());
	
	// Simplices incident to vertex
	Simplex * const *iSimplex = &incidentSimplices[iIndex[vertex->getLocalIndex()]];
	const int nIncident = iIndex[vertex->getLocalIndex()+1]-iIndex[vertex->getLocalIndex()];
	/*
	for (int i=0;i<nIncident;++i)
	  if (!iSimplex[i]->isLocal()) exit(-1);
	*/
	if (voxel==NULL)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>
	      ("While projecting: the mesh is not fully contained within the AMR grid.\n");
	    //throw std::runtime_error
	    //  ("Projecting: mesh is not fully contained within the AMR grid.");
	    exit(-1);
	  }
	
	// First, for each vertex compute its contribution to the voxel and retrieve 
	// the incident edges. For each vertex, only the subset of the edges
	// it owns is actually returned, other edges belonging to the vertex at 
	// their other extremity
	
	Base::getVertexContribAndEdges
	  (iSimplex,nIncident,vertex,weightFunctor,ownedEdges);

	// Add the current vertex contrib from the segments this vertex owns. 
	// The other extremity's contrib will be added when raytracing the segment
	//(*out)=std::make_pair(voxel,vertexContrib);++out;
	//++nProcessed;	

	// Now we can compute the contribution of the owned simplices edges intersection
	// with the voxels facets by raytracing them in the AMR grid	
	for (int i=0;i<ownedEdges.size();++i)
	  {
	    IncidentEdge &edge=ownedEdges[i];
	    // Note that destVoxel may be NULL, this should still be OK ...
	    Voxel *destVoxel = amr->getVoxelAt(edge.otherVertex->getCoordsConstPtr());

	    // First add the contribution of the current vertex	    
	    IFloat contribCoords[NDIM];
	    std::copy_n(edge.vertex->getCoordsConstPtr(),NDIM,contribCoords);
	    
	    if ((Pass==0)||(contribSumInterface.needReprojection(voxel)))
	      nProcessed+=Base::template addContrib<+1>
		(contribCoords,weightFunctor,edge,voxel,out);

#ifdef DEBUG_CHECK
	    if ((voxel->getIndex() == DEBUG_CHECK)&&
		(glb::mpiComWorld->rank()==DEBUG_CHECK_RANK)) 
	      {
#pragma omp critical
		glb::console->print<LOG_STD_ALL>("ADDING simplex vertex");// : %g\n",
		//hlp::numericStaticCast<double>(contrib));
	      }
#endif
	   
	    // and add the contribution of the other extremity of the segment
	    // if necessary
	    if (destVoxel!=NULL)
	      {		
		IFloat contribCoords[NDIM];
		std::copy_n(edge.otherVertex->getCoordsConstPtr(),NDIM,contribCoords);

		if ((Pass==0)||(contribSumInterface.needReprojection(destVoxel)))
		  nProcessed+=Base::template addContrib<-1>
		    (contribCoords,weightFunctor,edge,destVoxel,out);
	
#ifdef DEBUG_CHECK
		if ((destVoxel->getIndex()==DEBUG_CHECK)&&
		    (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK)) 
		  {
#pragma omp critical
		    glb::console->print<LOG_STD_ALL>
		      ("ADDING opposite simplex vertex");// : %g\n",
		      //hlp::numericStaticCast<double>(-contrib));
		  }	
#endif	
	      }

	    // now we can raytrace only if the segment is not entirely contained 
	    // within a single voxel
	    if (destVoxel == voxel) continue;
	    
	    raytracer.reset(edge.vertex->getCoordsConstPtr(),
			    edge.otherVertex->getCoordsConstPtr(),
			    voxel);
	    /*
	    bool dbg=false;
	    if ((destVoxel->getIndex()==162129587667468288)&&
		(voxel->getIndex()    ==162129587633913856))
	      {
		dbg=true;
		std::cout.precision(20);
		std::cout<<std::scientific;
	      }
	    if (dbg)
	      std::cout << "Raytracing from [" 
			<< edge.vertex->getCoordsConstPtr()[0] <<"," 
			<< edge.vertex->getCoordsConstPtr()[1] << "] to [" 
			<< edge.otherVertex->getCoordsConstPtr()[0] << "," 
			<< edge.otherVertex->getCoordsConstPtr()[1] << "]"<<std::endl;
	    if (dbg) std::cout << "Starting at voxel " << voxel->getIndex() << std::endl;
	    */
	    
	    Voxel *nextVoxel=NULL;	    
	    do
	      {	
		nextVoxel = raytracer.advance();
		/*
		if (dbg) {
		  const RayCoord *exitPoint=raytracer.getExitPointConstPtr();
		  std::cout << "Entering voxel " 
			    << nextVoxel->getIndex() << " at ["
			    << exitPoint[0]<<","<<exitPoint[1]<<"]"<<std::endl;
		}
		if (dbg) std::cout 
			   << "Normal along "<<raytracer.getExitFaceNormalDim()
			   << " with sign " << raytracer.getExitFaceNormalSign()<<std::endl;
		*/
		
#ifdef DEBUG_CHECK
		bool debugVoxelFound=false;
		
		if ((raytracer.getCurVoxel()->getIndex()==DEBUG_CHECK)&&
		    (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))
		  {
		    const RayCoord *exitPoint=raytracer.getExitPointConstPtr();
#pragma omp critical
		    {
		      printf("ADDING RAY CONTRIB (CURR) (normal: d=%d s=%d)\n",
			     hlp::numericStaticCast<int>(raytracer.getExitFaceNormalDim()),
			     hlp::numericStaticCast<int>(raytracer.getExitFaceNormalSign())
			     );
		      std::cout<<"exit point: ("<<
			exitPoint[0]<<" "<<exitPoint[1]<<" "<<exitPoint[2]<<")\n";
		    }
		    debugVoxelFound=true;
		  }
		if ((raytracer.getNextVoxel()->getIndex()==DEBUG_CHECK)&&
		    (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))
		  {
		    const RayCoord *exitPoint=raytracer.getExitPointConstPtr();
#pragma omp critical
		    {
		      printf("ADDING RAY CONTRIB (NEXT) (normal: d=%d s=%d)\n",
			     hlp::numericStaticCast<int>(raytracer.getExitFaceNormalDim()),
			     hlp::numericStaticCast<int>(raytracer.getExitFaceNormalSign())
			     );
		      std::cout<<"exit point: ("<<
			exitPoint[0]<<" "<<exitPoint[1]<<" "<<exitPoint[2]<<")\n";
		    }
		    debugVoxelFound=true;
		  }

		if (debugVoxelFound)
		  {
		    const Coord *start=edge.vertex->getCoordsConstPtr();
		    const Coord *stop=edge.otherVertex->getCoordsConstPtr();
		    printf("Ray is going from V%d(%20.20e,%20.20e,%20.20e) to V%d(%20.20e,%20.20e,%20.20e) (%d %d %d)(%d %d %d)\n",
			   edge.vertex->getLocalIndex(),start[0],start[1],start[2],
			   edge.otherVertex->getLocalIndex(),stop[0],stop[1],stop[2],
			   start[0]==int(start[0]),start[1]==int(start[1]),
			   start[2]==int(start[2]),stop[0]==int(stop[0]),
			   stop[1]==int(stop[1]),stop[2]==int(stop[2]));
		  }
#endif
		
		if ((Pass==0)||
		    (contribSumInterface.needReprojection(raytracer.getCurVoxel()))||
		    (contribSumInterface.needReprojection(raytracer.getNextVoxel())))
		  nProcessed += Base::addVoxelFacetSimplexEdgeContrib
		    (raytracer,weightFunctor,edge,out);
	
#ifdef DEBUGDUMP	
		const RayCoord *exitPoint=raytracer.getExitPointConstPtr();
		fprintf(fl[1],"%e %e %e\n",exitPoint[0],exitPoint[1],exitPoint[2]);
#endif
	
	      } while (nextVoxel != destVoxel);	     	
	  }
  
	// And finaly compute the voxels contributions (i.e. of their vertices and in 3D 
	// of the intersections of their segments with owned simplices facets)
	for (int i=0;i<nIncident;++i)
	  {
	    Simplex *simplex=iSimplex[i];

	    bool ownSimplex=true;
	    // A simplex is owned by a vertex if the vertex's address is the highest among 
	    // that of the local vertices in the simplex
	    for (int j=0;j<Simplex::NVERT;++j)
	      {
		Vertex *sv=simplex->getVertex(j);
		if ((vertex < sv)&&(sv->isLocal()))
		  ownSimplex=false;
	      }
	    
	    // we own this simplex, so we have to compute the contribution of the voxels 
	    // it intersects
	    if (ownSimplex) 
	      nProcessed+=addVoxelContribs<Base,Pass>(simplex,overlap,weightFunctor,out);
	  }
      }
    
    return nProcessed;
  }

  // This function computes the contributions from the intersection of a given simplex with 
  // the voxel's vertices as well as that from the voxel edges intersection with 
  // the simplex facets (the later is only needed in 3D)
  // TODO?: actually vertices on boundary of Union of overlapping voxels do not have to be 
  // tested as they cannot possibly belong to the simplices ! For that to work, we have
  // to make sure we are adding epsilon to each boundary of the bbox in getOverlap
  template <class Base, int Pass, class WF, class OutputIterator>
  long addVoxelContribs(Simplex *simplex, std::vector<Voxel*> &overlap, 
			const WF &weightFunctor, OutputIterator out) 
  {
    typedef typename Base::Float  IFloat;
    typedef typename Base::HFloat HFloat;
    typedef typename Simplex::Coord Coord;

    long nProcessed=0;  
    Coord simplexBBox[2][NDIM];       
    const Coord *refCoords=simplex->getVertex(0)->getCoordsConstPtr(); 

    FacetHandle ownedFacets[Simplex::NFACET];
    int ownedFacetsIndex[Simplex::NFACET];
    int ownedFacetsCount=0;

    if (NDIM!=2) 
      ownedFacetsCount=getOwnedFacets(simplex,ownedFacets,ownedFacetsIndex);

    if (weightFunctor.getTag(simplex)!=0)
      {
	int tag=weightFunctor.getTag(simplex);
	if (tag&ProjectionTag::skip)//(1<<2))
	  {
	    //skip this simplex
	  }
	else if (tag&ProjectionTag::sampleFromVertices)
	  {
	    double mass = weightFunctor.getMass(simplex)/(dimFactor*Simplex::NVERT);
	    for (int i=0;i<Simplex::NVERT;++i)
	      {
		const Coord *coords=simplex->getVertex(i)->getCoordsConstPtr(); 
		Voxel *v = amr->getVoxelAt(coords);		
		(*out)=std::make_pair(v,mass);++out;	
	      }
	    nProcessed+=Simplex::NVERT;
	  }
	else if (tag&ProjectionTag::sampleNoOverlap)//(1<<0))
	  {
	    
	    // This simplex is fully contained within a voxel
	    Voxel *v = amr->getVoxelAt(refCoords);
	    double mass = weightFunctor.getMass(simplex)/dimFactor;
	    (*out)=std::make_pair(v,mass);++out;
	    ++nProcessed;
	  }
	else if (tag&ProjectionTag::sampleWithOverlap)//(1<<1))
	  {	    
	    // This simplex spans 2+ voxels and needs to be sampled
	    static const long nGen=30;
	    static const long contribCapacity=nGen*100;
	    Coord samples[nGen][NDIM];
	    std::pair<Voxel*,double> contrib[contribCapacity];
	    std::vector< std::pair<Voxel*,double> > extraContrib;
	    std::unordered_set<ICoord> uniqueVoxels;	    	    
	    
	    long count=0;
	    do{
	      double sampleValues[nGen];
	      weightFunctor.sample(simplex,nGen,&samples[0][0],sampleValues);	
	      for (long i=0;i<nGen;++i)
		{		
		  Voxel *v = amr->getVoxelAt_noCheck(samples[i]);
		  
		  if (count < contribCapacity)
		    contrib[count]=std::make_pair(v,sampleValues[i]);
		  else 
		    extraContrib.push_back(std::make_pair(v,sampleValues[i]));
		  
		  ++count;		  
		  uniqueVoxels.insert(v->getIndex());
		}
	      //count+=nGen;
	    } while (count <= uniqueVoxels.size()*samplesPerVoxel);
	    
	    nProcessed += count;
	    const double factor = double(nGen)/(dimFactor*count);
	    for (long i=0;i<count-extraContrib.size();++i)
	      {
		contrib[i].second *= factor;
		(*out)=contrib[i];
		++out;
	      }
	    for (long i=0;i<extraContrib.size();++i)
	      {
		extraContrib[i].second *= factor;
		(*out)=extraContrib[i];
		++out;
	      }	  
	  }

	// If NDIM is 2, there is nothing left to do, but if NDIM==3, we still
	// have to check facets ...
	if (NDIM==2) return nProcessed;
	else 
	  {
	    int j=0;
	    // Remove null weighted facets
	    for (int i=0;i<ownedFacetsCount;++i)
	      {
		Simplex *nei=ownedFacets[i]->getOppositeSimplex();
		if ((nei==NULL)||
		    (!nei->isLocal())||
		    (weightFunctor.getTag(nei)==0))	
		  {
		    if (i!=j) 
		      {
			ownedFacets[j]=ownedFacets[i];
			ownedFacetsIndex[j]=ownedFacetsIndex[i];
		      }
		    j++;
		  }		    		  
	      }
	    ownedFacetsCount=j;
	    if (ownedFacetsCount==0) return nProcessed;
	  }
      }
          
    // First compute the bounding box of the simplex    
    for (int j=0;j<NDIM;++j) 
      simplexBBox[0][j]=simplexBBox[1][j]=refCoords[j];
    
    bool simplexCoordsAreConsistent=true;

    for (int j=1;j<Simplex::NVERT;++j)
      {
	Vertex *v=simplex->getVertex(j);
	const Coord *coords=v->getCoordsConstPtr();
	for (int k=0;k<NDIM;++k)
	  {
	    // We have to be a bit carefull with periodic boundaries !
	    Coord c = meshGeometry->checkCoordConsistency(coords[k],refCoords[k],k);
	    
	    if (simplexBBox[0][k]>c)
	      simplexBBox[0][k]=c;
	    if (simplexBBox[1][k]<c)
	      simplexBBox[1][k]=c;

	    if (c!=coords[k]) 
	      {		
		simplexCoordsAreConsistent=false;		
	      }
	  }	
      }
    
    overlap.clear();
    amr->getLeavesBBoxOverlap(simplexBBox,std::back_inserter(overlap));

    // This is used only to compute the voxel segments and 
    // simplex facets intersections in 3D
    IFloat facetNormal[Simplex::NFACET][NDIM];
    if (NDIM>2)
      {	
	// Compute normalized facet normals 
	for (int f=0;f<ownedFacetsCount;++f)
	  {
	    ownedFacets[f]->template computeProjectedNormal<IFloat*,MGP,IFloat>
	      (&facetNormal[f][0],meshGeometry);
	    meshGeometry->template normalize_noCheck<IFloat,NDIM>(facetNormal[f]);	    
	  }
	/*
	// Compute normalized facet normals 
	for (int f=0;f<Simplex::NFACET;++f)
	  {
	    Simplex *nei=simplex->getNeighbor(f);
	    if ((nei < simplex)||(!nei->isLocal()))
	      {
		simplex->getFacetHandle(f)->
		  computeProjectedNormal(facetNormal[f],meshGeometry);
		meshGeometry->template normalize_noCheck<Coord,NDIM>(facetNormal[f]);
	      }
	  }
	*/
      }

    Coord corner[2][NDIM];
    Coord coords[NDIM];
    Coord coords2[NDIM];

    Voxel *neighbors[Voxel::NNEI];
    std::pair<int,int> segments[Voxel::NSEG];
    //std::vector< std::pair<int,int> > segments;

    for (unsigned long ovi=0;ovi<overlap.size();ovi++)
      {
	Voxel *curOverlap=overlap[ovi];
	// get the voxel corners coordinates
	amr->index2CornerCoordsAndOpp(curOverlap->getIndex(),
				      curOverlap->getLevel(),
				      corner[0],corner[1]);

	// We only compute voxels' segments intersection with simplices facets in 3D
	if (NDIM>2)
	  {	
	    int nSegs = amr->getOwnedSegments(curOverlap,segments);
	    for (int s=0;s<nSegs;++s)
	      {
		int vertexId = segments[s].first;		
		int dim=segments[s].second;

		if (Pass>0)
		  {
		    amr->getSegmentNeighbors(curOverlap,vertexId,dim,neighbors);
		    bool dontcare=true;
		    for (int i=0;i<Voxel::NSEGNEI;++i)
		      if (contribSumInterface.needReprojection(neighbors[i]))
			dontcare=false;
		    if (dontcare) continue;
		  }
		
		// Compute the coordinates of the vertex from the voxel corners	
		for (int d=0;d<NDIM;d++)
		  {
		    coords2[d]=corner[(vertexId>>d)&1][d];
		    coords[d]=coords2[d];
		  }
	
		// We precompute whether coordinates are consistent to avoid having to 
		// check it every time we test for intersection
		bool coordsAreConsistent=simplexCoordsAreConsistent;

		if (coords[dim] != meshGeometry->
		    checkCoordConsistency(coords[dim],refCoords[dim],dim))
		  coordsAreConsistent=false;
		
		// Check whether the segment is entirely outside the bbox as in that
		// case we dont need to test it for intersection (which is more expensive)
		if (dim==0)
		  {
		    Coord c1 = meshGeometry->
		      checkCoordConsistency(coords[1],refCoords[1],1);
		    Coord c2 = meshGeometry->
		      checkCoordConsistency(coords[2],refCoords[2],2);
		    if ((c1<simplexBBox[0][1])||
			(c2<simplexBBox[0][2])||
			(c1>simplexBBox[1][1])||			
			(c2>simplexBBox[1][2]))
		      continue;
		    else if ((c1!=coords[1])||(c2!=coords[2]))
		      coordsAreConsistent=false;
		  }
		else if (dim==1)
		  {
		    Coord c1 = meshGeometry->
		      checkCoordConsistency(coords[0],refCoords[0],0);
		    Coord c2 = meshGeometry->
		      checkCoordConsistency(coords[2],refCoords[2],2);
		    if ((c1<simplexBBox[0][0])||
			(c2<simplexBBox[0][2])||
			(c1>simplexBBox[1][0])||			
			(c2>simplexBBox[1][2]))
		      continue;
		    else if ((c1!=coords[0])||(c2!=coords[2]))
		      coordsAreConsistent=false;
		  }
		else
		  {
		    Coord c1 = meshGeometry->
		      checkCoordConsistency(coords[0],refCoords[0],0);
		    Coord c2 = meshGeometry->
		      checkCoordConsistency(coords[1],refCoords[1],1);
		    if ((c1<simplexBBox[0][0])||
			(c2<simplexBBox[0][1])||
			(c1>simplexBBox[1][0])||			
			(c2>simplexBBox[1][1]))
		      continue;
		    else if ((c1!=coords[0])||(c2!=coords[1]))
		      coordsAreConsistent=false;
		  }
			
		coords2[dim]=corner[1][dim];		

		// If we are not periodic, we consider that a facet lying on the boundary
		// of the bounding box is actually inside it.
		// Note that the two tests are skipped for periodic boundaries.
		if (amrGeometry->onBoundary(coords2)) 
		  {
		    for (int i=0;i<NDIM;++i)
		      {
			if (coords2[i]<=amr->getBBoxMin(i)) coords2[i]-=amr->getEpsilon(i);
			if (coords2[i]>=amr->getBBoxMax(i)) coords2[i]+=amr->getEpsilon(i);
		      }		    
		  }
		if (amrGeometry->onBoundary(coords)) 
		  {
		    for (int i=0;i<NDIM;++i)
		      {
			if (coords[i]<=amr->getBBoxMin(i)) coords[i]-=amr->getEpsilon(i);
			if (coords[i]>=amr->getBBoxMax(i)) coords[i]+=amr->getEpsilon(i);
		      }		    
		  }

		bool alreadyCrossed=false;
		IFloat contribSign[NDIM];
		//for (int f=0;f<Simplex::NFACET;++f)
		for (int f=0;f<ownedFacetsCount;++f)
		  {	
		    /*
		    if (ownedFacets[f]->getSimplex()->getLocalIndex()==1414079)
		      {
			glb::console->print<LOG_STD_ALL>
			  ("TESTING voxel edge : @D%d(%.10lg %.10lg %.10lg)\n",
			   dim,coords[0],coords[1],coords[2]);
			ownedFacets[f]->template print<LOG_STD>();
		      }
		    */
		    /*
		    Simplex *nei=simplex->getNeighbor(f);
		    if ((nei>simplex)&&(nei->isLocal())) continue;
		    */
		    //if (!simplexOwnsFacet(simplex,f)) continue;

		    IFloat crossCoord;			    
		    if (ownedFacets[f]->template intersectSegment<MGP,IFloat,IFloat>
			(coords,dim,coords2[dim],meshGeometry,
			 crossCoord,coordsAreConsistent))
		      {					
			if (!alreadyCrossed)
			  {
			    amr->getSegmentNeighbors(curOverlap,vertexId,dim,neighbors);
			    for (int i=0;i<NDIM;++i)
			      contribSign[i]=segmentNeighborDir[vertexId][dim][i];
			    //amr->getSegmentNeighborsDirection(vertexId,dim,contribSign);

			    alreadyCrossed=true;
			  }
			
			// Compute the coordinates of the intersection point
			IFloat intersection[NDIM];
			for (int i=0;i<NDIM;++i) 
			  intersection[i]=coords[i];
			intersection[dim]=crossCoord;

			/* * DEBUG * */
			/*
			int dbg=0;
			for (int i=0;i<NDIM;++i)
			if ((ownedFacets[f]->getVertex(i)->getLocalIndex()==2592)||
			    (ownedFacets[f]->getVertex(i)->getLocalIndex()==260575)||
			    (ownedFacets[f]->getVertex(i)->getLocalIndex()==260576))
			  dbg++;
			if (dbg!=NDIM) dbg=0;
			if (dbg)
			  {
			    std::cout<<"Facet intersection FOUND @"<< dim <<"("
				     <<intersection[0]<<" "
				     <<intersection[1]<<" "
				     <<intersection[2]<<")"<<std::endl;
			    ownedFacets[f]->template print<LOG_STD>();
			  }
			*/
			/* * END * */
#ifdef DEBUGDUMP	    
			fprintf(fl[0],"%e %e %e\n",intersection[0],intersection[1],intersection[2]);
#endif

#ifdef DEBUG_CHECK
			for (int i=0;i<Voxel::NSEGNEI;++i)
			  {
			    if ((neighbors[i]->getIndex() == DEBUG_CHECK)&&
				(glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))
			      {
#pragma omp critical
				{
				  std::cout << "ADDING voxel edge : @("
					    << intersection[0] <<" "
					    << intersection[1] <<" "
					    << intersection[2] <<")D="<<dim
					    <<"; sign=["
					    <<contribSign[0]<<","
					    <<contribSign[1]<<","
					    <<contribSign[2]<<"]"<<std::endl;
				  std::cout << "Normal: (" <<facetNormal[f][0]<<","
					    <<facetNormal[f][1]<<","
					    <<facetNormal[f][2]<<")"<<std::endl;
				
				  /*
				    glb::console->print<LOG_STD_ALL>
				    ("ADDING voxel edge : @(%.10lg %.10lg %.10lg)D=%d\n",
				    intersection[0],intersection[1],intersection[2],dim);
				  */
				  ownedFacets[f]->template print<LOG_STD>();
				}
				  
			      }
			  }
#endif

			bool dbg=false;
			/*
			if ((ownedFacets[f]->getSimplex()->getLocalIndex()==549310)&&
			    (ownedFacets[f]->getOppositeVertexIndex()==1))
			  dbg=true;
			*/

			nProcessed += Base::addVoxelEdgeSimplexFacetContrib
			  (simplex,ownedFacetsIndex[f],facetNormal[f],dim,
			   contribSign,intersection,neighbors,
			   weightFunctor,out,coordsAreConsistent,dbg);		
		      }//if (intersect)
		  }//foreach simplex facets
	      }//foreach voxel segments
	  }//if (NDIM>2)
	
	// If the simplex is tagged, its mass weight and grad is 0 so we can return
	if (weightFunctor.getTag(simplex)!=0) continue;

	// Now test if voxel vertices are inside the simplex
	for (int vertexId=0;vertexId<AMR::Voxel::NVERT;++vertexId) 
	  {	    
	    // and proceed only if the voxel owns this vertex ...
	    if (curOverlap->getVertexFlag(vertexId)) 
	      {	
		bool onOverlapBoundary=false;  
		bool coordsAreConsistent=simplexCoordsAreConsistent;
		// Compute the coordinates of the vertex from the voxel corners
		// FIXME: this depends on localAMR conventions, should be in localAmrGrid ?
		for (int d=0;d<NDIM;d++)
		  {
		    coords[d]=corner[(vertexId>>d)&1][d];
		    Coord c = meshGeometry->checkCoordConsistency(coords[d],refCoords[d],d);
		    if ((c<simplexBBox[0][d]-amr->getEpsilon(d))||
			(c>simplexBBox[1][d]+amr->getEpsilon(d)))
		      onOverlapBoundary=true;
		    else if (c!=coords[d])
		      coordsAreConsistent=false;
		  }
		if (onOverlapBoundary) continue;
		
		if (simplex->pointIsInside(coords,meshGeometry,coordsAreConsistent))
		  {
		    amr->getVertexNeighbors(curOverlap,vertexId,neighbors);
		    
		    // Always compute the contribution of the voxels vertices.
		    // NOTE: because of simulation of simplicity, voxel's vertices on 
		    // the left boundary cannot fall inside any simplex, so their contrib 
		    // is 0, but those on the right boundary may. 
		    // We previously wrongly counted the contrib of simplex vertices that
		    // fall exactly on a right boundary, so we compensate for that by 
		    // preventing the voxel's boundary vertices from contributing. 
		    // For periodic boundary, this does not change anything as there is no
		    // boundary anyway !
		    if (!amrGeometry->onBoundary(coords))
		      {	   
#ifdef DEBUG_CHECK
			for (int i=0;i<Voxel::NNEI;++i)
			  {
			    if ((neighbors[i]->getIndex() == DEBUG_CHECK)&&
				(glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))
			      {		
#pragma omp critical
				{
				  glb::console->print<LOG_STD>("ADDING VERTEX CONTRIB (%d) @ %20.20e %20.20e %20.20e\n",vertexId,coords[0],coords[1],coords[2]);
				//glb::console->print<LOG_STD>("Local:%d\n",simplex->isLocal());
				  simplex->template print<LOG_STD>();
				}
			      }
			  }
#endif
			IFloat contribCoords[NDIM];
			std::copy_n(coords,NDIM,contribCoords);
			nProcessed += Base::
			  addVoxelCornerContrib(contribCoords,vertexNeighborDir[vertexId],
						neighbors,weightFunctor,
						simplex,out,coordsAreConsistent);
		      } // onBoundary
		  } // isInside
	      } // vertexFlags
	  } // vertexId
	  
      } // overlap
    
    return nProcessed;
  }
  
private:
  ContribSumInterface contribSumInterface;
  //HPBase highPrecisionBase;

  int getOwnedFacets(Simplex *simplex, FacetHandle *fh, int *index)
  {
    // We only need to check each facet once, not once for each simplex 
    // it is adjacent to. The convention is that it belongs to the adjacent 
    // local simplex with highest adress ... These also ensure that we do not
    // forget boundary facets ...

    int nOwned=0;
    for (int i=0;i<Simplex::NNEI;++i)
      {
	Simplex *nei=simplex->getNeighbor(i);
	if ((nei<simplex)||(!nei->isLocal()))
	  {
	    fh[nOwned]=simplex->getFacetHandle(i);
	    index[nOwned]=i;
	    ++nOwned;	    
	  }
      }
    return nOwned;
  }
  
  bool simplexOwnsFacet(Simplex *simplex, int index)
  {
    Simplex *nei=simplex->getNeighbor(index);
    return ((nei<simplex)||(!nei->isLocal()));
  }

  // Acts like a std::back_insert( container< std::pair<T*,T> > ) but instead of pushing 
  // the std::pair<T*,T> p in a container when assigned to, do (*p.first)+=p.second
  template <class CSI>
  class FakeBackInserterT
  {
    typedef FakeBackInserterT<CSI> MyType;
  public:
    FakeBackInserterT(CSI &csi_):csi(csi_)
    {}
    
    template <class VP>
    MyType& operator=(const VP &rhs)
    {       
      //if (rhs.second != rhs.second) throw 0;
      //rhs.first->data += rhs.second;
      csi.add(rhs);
     
#ifdef DEBUG_CHECK
      static int count=0;
      long vid=DEBUG_CHECK;
      //long vid=4602678820611293184L;
      
      if ((rhs.first->index==vid)&&
	  (glb::mpiComWorld->rank()==DEBUG_CHECK_RANK))//&&(rhs.second!=0))
	{
	  //const typename CSI::SumType &val=rhs.second;
	  //double res = val.template convert_to<double>();
	  /*
	  double res = hlp::numericStaticCast<double>(rhs.second);//static_cast<double>(val);
	  double res2= hlp::numericStaticCast<double>(csi.getValue(rhs));
	  printf("%2.2d: Adding %16.16e to voxels %16.16ld -> %16.16e\n",++count,res,vid,res2);	  
	  */
#pragma omp critical
	  {
	    std::cout.precision(20);
	    std::cout<<std::scientific<<++count<<": Adding "<<rhs.second
		     <<" to voxels "<<vid
		     <<" => "<<csi.getValue(rhs)<<std::endl;
	  }
	}    
#endif

      return *this; 
    }

    MyType& operator *() { return *this; }
    MyType& operator ++() { return *this; }  
  private:
    CSI &csi;
  };
  
  template <bool TO_DENSITY>
  class CorrectDimFactorVisitorT
  {
    typedef typename AMR::Voxel Voxel;    
    typedef typename hlp::IsTrueT<TO_DENSITY>::Result ToDensity;
    typedef typename AMR::Data Data;
  public:
    CorrectDimFactorVisitorT(AMR *amr_, Data factor_):amr(amr_),factor(factor_)
    {}

    static void initialize(Voxel *rootVoxel) {}

    bool visit(Voxel *voxel) const
    {
      if (voxel->isLeaf())
	{	 
	  // if (voxel->data != 0)
	  //   voxel->print(amr,"DATA:");	  
	  
	  voxel->data *= getFactor(voxel,ToDensity());
	  return false;
	}   
      return true;
    }    

    static bool visit(Voxel *voxel,int i) 
    {
      return true;
    }

    static void visited(const Voxel *voxel) {}

    static void visited(Voxel *voxel, int i) 
    {
    }
    
  private:
    Data getFactor(Voxel *voxel, hlp::IsTrue) const
    {return amr->getVoxelInverseVolume(voxel->getLevel())*factor;}
    Data getFactor(Voxel *voxel, hlp::IsFalse) const
    {return factor;}

    const AMR *amr;
    const Data factor;
  };
  
private:
  // local helper functions for debugging ...
  template <class F>
  bool findFacet(double x1,double y1,double z1,
		 double x2,double y2,double z2,
		 double x3,double y3,double z3,
		 const F &facet, double tol=2.E-5)
  {
    double coord[3][3]={{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}};

    int count=0;   
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	{	  
	  int k=0;
	  for (;k<3;++k)
	    if (fabs(facet->getVertex(i)->getCoord(k)-coord[j][k])>tol)
		break;
	  if (k==3) {count++;break;}
	}

    return (count==3);
  }

  bool findVertex(double x1,double y1,double z1,double vx, double vy, double vz, double tol=2.E-5)
  {
    if ((fabs(meshGeometry->correctCoordsDiff(vx-x1,0))<=tol)&&
	(fabs(meshGeometry->correctCoordsDiff(vy-y1,1))<=tol)&&
	(fabs(meshGeometry->correctCoordsDiff(vz-z1,2))<=tol))
      {
	return true;
      }
    return false;
  }
  //long countT;
};

/** \}*/

#ifdef DEBUGDUMP
#undef PRINTCOORD
#undef PRINTVERTEX
#undef DEBUGDUMP
#endif

#ifdef DEBUG_CHECK
#undef DEBUG_CHECK_RANK
#undef DEBUG_CHECK_COORD
#undef DEBUG_CHECK
#endif

#include "../internal/namespace.footer"

#endif
