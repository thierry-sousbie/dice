#ifndef __COLDICE_POST_HXX__
#define __COLDICE_POST_HXX__

#include "../coldice.hxx"
// This is just a convient hack to compute statistics from a snapshot. 

template <class VlasovPoissonSolver>
class ColdicePost : public VlasovPoissonSolver
{
public:
  typedef VlasovPoissonSolver Base;
  typedef ColdicePost<VlasovPoissonSolver> MyType;
 
  typedef typename Base::Mesh Mesh;
  typedef typename Base::MeshParams MeshParams;
  typedef typename Base::Simplex Simplex;

  typedef typename Base::SimplexFunctor SimplexFunctor;
  typedef typename Base::VertexFunctor  VertexFunctor;
 
  static const int NDIM          = Base::NDIM;
  static const int NDIM_W        = Base::NDIM_W;
  static const int BOUNDARY_TYPE = Base::BOUNDARY_TYPE;

  static std::string parserCategory() {return "post";}
 
  template <class SP, class R, class PM>
  ColdicePost(MeshParams &meshParams,
	      SP &solverInterfaceParams,
	      R *reader,
	      PM& paramsManager,
	      float serializedVersion,
	      dice::MpiCommunication *mpiCom_):
    Base(meshParams,solverInterfaceParams,reader,
	 paramsManager,serializedVersion,mpiCom_)
  {     
    typename PM::Parser *parser = paramsManager.getParser();
     
    gridLevel=Base::fftGridLevel;
    gridLevel=parser->
      get("gridLevel",parserCategory(),gridLevel,
	  "Sets the resolution of the density grid");

    projectionOrder=Base::projectionOrder;
    projectionOrder=parser->
      get("projectionOrder",parserCategory(),projectionOrder,
	  "Order of the exact mass projection (0 or 1)");

    refineMesh=0;
    refineMesh=parser->
      get("refineMesh",parserCategory(),refineMesh,
	  "Refine the mesh using tracers.");

    dumpAmr = 0;
    dumpAmr = parser->
      get("dumpAmr",parserCategory(),dumpAmr,
	  "Whether to dump AMR files.");

    dumpDensity = 0;
    dumpDensity = parser->
      get("dumpDensity",parserCategory(),dumpDensity,
	  "Whether to dump density grid files");

    dumpMesh = 0;    
    dumpMesh = parser->
      get("dumpMesh",parserCategory(),dumpMesh,
	  "Whether to dump the unstructured mesh files");   

    dumpSubmesh = 0;    
    dumpSubmesh = parser->
      get("dumpSubmesh",parserCategory(),dumpSubmesh,
	  "Whether to dump the unstructured mesh subsets files");

    dumpRadialGridDensity = 0;  
    dumpRadialGridDensity = parser->
      get("dumpRadialGridDensity",parserCategory(),dumpRadialGridDensity,
	  "Whether to dump the radial grid density (from projected grid)");

    dumpRadialMeshDensity = 0;  
    dumpRadialMeshDensity = parser->
      get("dumpRadialMeshDensity",parserCategory(),dumpRadialMeshDensity,
	  "Whether to dump the radial mesh density (from the mesh)");

    radialMeshDensitySamples = 10;  
    radialMeshDensitySamples = parser->
      get("radialMeshDensitySamples",parserCategory(),radialMeshDensitySamples,
	  "The number of samples to use per overlap for radial mesh density computation");
    
    
    needProjectedDensity=(dumpDensity||dumpAmr||dumpRadialMeshDensity);  
  }
  
  template <class SP>
  void onInitialize(Mesh *m, const SP& solverParams, double t, bool restart)
  { 
    if (!(dumpMesh||dumpDensity||dumpAmr||dumpRadialGridDensity||dumpRadialMeshDensity))
      {
	dice::glb::console->printFlush<dice::LOG_STD>("Post : nothing to do ...\n");
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("At least one of the following options must be set :\n");
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("    -post.dumpMesh -post.dumpAmr -post.dumpDensity.\n");

	exit(0);
      }

    Base::init_impl(m,solverParams,t,restart,(needProjectedDensity)?gridLevel:-1);     
  }

   void onNewTimeStep(long stepId, double curTime, double &newDeltaT)
  { 
    Base::onNewTimeStep(stepId,curTime,newDeltaT);
    process();
    exit(0);
  }

  double checkRefine_getValue(Simplex *s)
  {
    Simplex *p=s->getPartner();
    if ((p!=NULL)&&((s->getLocalIndex() >= unrefinedSimplicesCount)||
		    (p->getLocalIndex() >= unrefinedSimplicesCount)))
      return 0;    
    
    return
      dice::slv::refine::poincareInvariantWithSegTracers_order1<Mesh>
      (s,Base::geometry).first;
  }
 
protected:
  
  void process()
  {
    double elapsed;

    dice::glb::console->printFlush<dice::LOG_STD>("Post treating :\n");
    dice::glb::console->indent();

    if (refineMesh)
      {
	dice::glb::console->printFlush<dice::LOG_STD>("Refining the mesh ... ");
	typename dice::TimerPool::Timer refineTimer;
	refineTimer.start();
	unrefinedSimplicesCount = Base::mesh->getNSimplices();
	Base::mesh->refine(this,dice::glb::num_omp_threads,0);
	elapsed=refineTimer.stop();
	dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed); 
      }

    dice::glb::console->printFlush<dice::LOG_STD>("Updating projected density ... ");    
    Base::updateDensityTimer->start();
    Base::updateProjectedDensityField();
    elapsed=Base::updateDensityTimer->stop();
    dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed); 

    if (dumpMesh||dumpSubmesh) Base::dumpNetwork(dumpMesh,dumpSubmesh);

    if (dumpRadialMeshDensity) Base::dumpMeshDensityProfile(true,radialMeshDensitySamples);
  
    if (needProjectedDensity)
      {
	Base::buildAmrGrid(gridLevel);
	/*
	dice::glb::console->printFlush<dice::LOG_STD>("Building adaptive grid from mesh (ML:%d) ... ",gridLevel);
	Base::amrBuildTimer->start();
	int level=gridLevel-Base::LocalAmrGrid::ROOT_LEVEL;

	if (Base::fastAmrBuild)
	  Base::localAmrDensity.buildFromMesh_Fast
	    (Base::mesh,1.0,level,true,dice::glb::num_omp_threads,true,NDIM>2);
	else
	  Base::localAmrDensity.buildFromMesh(Base::mesh,1.0,level,true);

	elapsed=Base::amrBuildTimer->stop();
	dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed); 
	*/

	// Project the local mesh's projectedDensity onto the AMR grid
	Base::projectTimer->start(); 
	if (Base::pForceHighResolution)
	  do_project<typename Base::HFloat,typename Base::HFloat,typename Base::HFloat,
		     false>();
	else
	  do_project<typename Base::IFloat,typename Base::HFloat,typename Base::SFloat>();

	elapsed=Base::projectTimer->stop();   

	if (dumpAmr)
	  Base::localAmrDensity.toVtk
	    (Base::fileDumps.getLocalFName(FileDumps::Amr).c_str());	
    
	if (dumpDensity||dumpRadialGridDensity)
	  {
	    Base::scatterDensityTimer->start();
	    dice::glb::console->printFlush<dice::LOG_STD>("Scattering density field (AMR->grid):\n");
	    // This is necessary when we convolve with poisson kernel inplace 
	    // and it would not hurt otherwise :)
	    Base::potential.reinterpretNFields(); 
	    Base::potential.setName("density");    
	    Base::potential.erase();
	    Base::potential.addAmrGrid(&(Base::localAmrDensity)); 
	    elapsed = Base::scatterDensityTimer->stop();
	    dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);
	    
	    if (dumpRadialGridDensity)
	      Base::dumpGridDensityProfile(true);		
	    
	    if (dumpDensity)
	      {
		Base::potential.toVtk
		  (Base::fileDumps.getGlobalFName(FileDumps::Density).c_str(),
		   Base::fileDumps.getLocalFNameFormat(FileDumps::Density).c_str());
	      }

	  }   
      }
    dice::glb::console->unIndent();
  } 
  
private:
  template <typename IT, 
	    typename HT, 
	    typename ST,
	    bool checkAccuracy=D_ENABLE_ACCURACY_CHECKING>
  void do_project()
  {
    if (projectionOrder==0)
      {
	const SimplexFunctor &projectedDensityFunctor=
	  *Base::mesh->getSimplexFunctorPtr("projectedDensity");  
	Base::localAmrDensity.template projectMesh0
	  <Mesh,SimplexFunctor,checkAccuracy,IT,HT,ST>
	  (Base::mesh,projectedDensityFunctor,
	   Base::accuracyLevel,
	   (Base::fastAmrBuild)?true:false, // Use simplex tags
	   dice::glb::num_omp_threads,
	   true);
      }
    else
      {
	const VertexFunctor &projectedDensityFunctor=
	  *Base::mesh->getVertexFunctorPtr("projectedDensity");
	const SimplexFunctor &projectedDensityGradientFunctor=
	  *Base::mesh->getSimplexFunctorPtr("projectedDensityGradient");  
	
	Base::localAmrDensity.template projectMesh1
	  <Mesh,VertexFunctor,SimplexFunctor,checkAccuracy,IT,HT,ST>
	  (Base::mesh,projectedDensityFunctor,projectedDensityGradientFunctor,
	   Base::accuracyLevel,
	   (Base::fastAmrBuild)?true:false, // Use simplex tags
	   dice::glb::num_omp_threads,
	   true); 
      }
  }

  int gridLevel;
  int projectionOrder;
  int refineMesh;
  int dumpMesh;
  int dumpSubmesh;
  int dumpAmr;
  int dumpDensity;

  int dumpRadialMeshDensity;
  int dumpRadialGridDensity;
  int radialMeshDensitySamples;

  double unrefinedSimplicesCount;

  bool needProjectedDensity;
};

#endif
