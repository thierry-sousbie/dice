#ifndef __COLDICE_HXX__
#define __COLDICE_HXX__

#include <limits> 
//#include <random>

#include "units.hxx"
#include "coldice_cellData.hxx"
#include "coldice_functors.hxx"
#include "coldice_NDnetFilters.hxx"
#include "coldice_meshQuadratureFunctors.hxx"
#include "coldice_gridQuadratureFunctors.hxx"
#include "coldice_iterators.hxx"
#include "coldice_fileDumps.hxx"
#include "cflCondition_type.hxx"

#include "init/initialConditions.hxx"

#include <dice/dice_globals.hxx>
#include <dice/tools/wrappers/boostMultiprecisionFloat128.hxx>
#include <dice/tools/wrappers/boostMultiprecisionQD.hxx>
#include <dice/mesh/mesh.hxx>
#include <dice/mesh/mesh_traits.hxx>
#include <dice/AMR/localAmrGrid.hxx>  
#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/IO/console.hxx>
#include <dice/tools/sort/ompPSort.hxx>
#include <dice/tools/sort/peanoHilbert.hxx>
#include <dice/tools/wrappers/likwidWrapper.hxx>
#include <dice/tools/helpers/objectList.hxx>
#include <dice/solver/solverTools.hxx>
#include <dice/grid/regularGrid.hxx>
#include <dice/FFT/FFTWConvolver.hxx>
#include <dice/geometry/barycentricCoordinates.hxx>
#include <dice/cosmo/poissonKernel.hxx>
//#include <cosmo/cosmology.hxx>

#include "config.h"

template <int D, int BT=dice::BoundaryType::PERIODIC>
class Coldice
{
public:
  typedef Coldice<D,BT> MyType;
  typedef dice::hlp::IsEnabled RefinementStatus; // enable refinement
  typedef dice::hlp::IsEnabled CoarseningStatus; // enable coarsening
  
  typedef dice::MeshTraitsT<D,2*D,BT,
			    VERTEX_REFINE_COORDS_METHOD,
			    ColdiceVertexData<2*D>,
			    ColdiceSimplexData<2*D> > MeshTraits;  
  
  typedef dice::MeshT<MeshTraits> Mesh;
  typedef dice::LocalAmrGridT<D,double,BT,D_AMR_ROOT_LEVEL> LocalAmrGrid;
  typedef dice::RegularGridT<D,double,BT> RegularGrid;
  typedef dice::FFTWConvolverT<RegularGrid> FFTSolver;

  static const int Periodic=(BT==dice::BoundaryType::PERIODIC);

  static const int POTENTIAL_INTERP_ORDER=2;

#ifdef D_FFT_POISSON_SOLVER_COMPUTE_DISPLACEMENT
  typedef dice::cosmo::PoissonKernelT<FFTSolver, dice::cosmo::pkDisplacement, Periodic> 
  FFTPoissonKernel;
  
  typedef dice::gridKernel::InterpolateT<D,1> 
  InterpolationKernel;
#else
  typedef dice::cosmo::PoissonKernelT<FFTSolver, dice::cosmo::pkPotential, Periodic> 
  FFTPoissonKernel;
    
  typedef dice::gridKernel::InterpolateCentralDiffT<D,POTENTIAL_INTERP_ORDER> 
  InterpolationKernel;
#endif

  // Used for projection 
  typedef D_PROJECTION_FLOAT_TYPE IFloat;
  typedef D_PROJECTION_HR_FLOAT_TYPE SFloat;
  typedef D_PROJECTION_HR_FLOAT_TYPE HFloat;

  // typedef dice::DoubleDouble IFloat;// Used for projection (for regular computations) 
  // typedef dice::DoubleDouble SFloat;// Used for projection (for summing contributions)
  // typedef dice::DoubleDouble HFloat;// Used for projection only (high precision)

  //typedef long double IFloat; // Used for projection (for regular computations)  
  //typedef long double SFloat; // Used for projection (for summing contributions)  
  //typedef dice::DoubleDouble HFloat; // Used for projection only (high precision)
  //typedef dice::Float128OrMore HFloat; // This is twice slower and slightly more precise
  //typedef dice::QuadDouble HFloat; // Extremely high precision !

  typedef InitialConditionsListT<Mesh> InitialConditionsList;
  typedef InitialConditionsFactoryT<InitialConditionsList> InitialConditionsFactory;
  typedef typename InitialConditionsFactory::IC IC;

  typedef typename RegularGrid::LocalGrid LocalGrid;
  typedef typename RegularGrid::Params RegularGridParams;

  typedef typename Mesh::Params        MeshParams;
  typedef typename Mesh::Simplex       Simplex;
  typedef typename Mesh::SegmentHandle SegmentHandle;  
  typedef typename Mesh::FacetHandle   FacetHandle;
  typedef typename Mesh::Vertex        Vertex;
  typedef typename Mesh::Coord         Coord;

  typedef typename Mesh::SimplexFunctor SimplexFunctor;
  typedef typename Mesh::VertexFunctor  VertexFunctor;

  typedef typename Mesh::GeometricProperties GeometricProperties;

  typedef typename Mesh::simplexPtr_iterator       simplexPtr_iterator;
  typedef typename Mesh::vertexPtr_iterator        vertexPtr_iterator;
  typedef typename Mesh::const_simplexPtr_iterator const_simplexPtr_iterator;
  typedef typename Mesh::const_vertexPtr_iterator  const_vertexPtr_iterator;
  
  typedef typename Mesh::simplexPtr_LG_iterator    simplexPtr_LG_iterator;
  typedef typename Mesh::vertexPtr_LG_iterator     vertexPtr_LG_iterator;
  typedef typename Mesh::simplexPtr_LGS_iterator   simplexPtr_LGS_iterator;
  typedef typename Mesh::vertexPtr_LGS_iterator    vertexPtr_LGS_iterator;

  typedef typename Mesh::CheckRefineReturnType     CheckRefineReturnType;

  typedef MeshAndTracersCoordsT<Mesh> MeshAndTracersCoords;
  typedef typename MeshAndTracersCoords::iterator MeshAndTracersCoords_iterator;

  typedef typename Simplex::Neighborhood SimplexNeighborhood;

  static const int NDIM          = Simplex::NDIM;
  static const int NDIM_W        = Simplex::NDIM_W;
  static const int BOUNDARY_TYPE = Mesh::BOUNDARY_TYPE;

  static std::string parserCategory() {return "solver";}
  static std::string classHeader() {return "vlasov_poisson_solver";}
  static float classVersion() {return 0.23;}
  static float compatibleSinceClassVersion() {return 0.17;}

  template <class SP, class R, class PM>
  Coldice(MeshParams &meshParams,
	  SP &solverInterfaceParams,
	  R *reader,
	  PM& paramsManager,
	  float serializedVersion,
	  dice::MpiCommunication *mpiCom_):
    mesh(NULL),
    geometry(NULL),  
    localAmrDensity(),
    clonedDensity("projectedDensity"),
    potential("potential"),
    mpiCom(mpiCom_),
    repartStatus(false)
  {       
    MyType::serializedVersion=serializedVersion;
   
    // First get the initial conditions type
    std::string icType=
      (serializedVersion<0.145)?
      std::string("uniformGrid"):
      std::string("sineWaves");

    std::stringstream icComment;
    icComment<< "The type of initial condition: " << InitialConditionsFactory::getList();
    icType=paramsManager.
      get("type","init",icType,reader,
	  PM::FILE_FIRST,
	  icComment.str(),
	  serializedVersion>0.145);
    ic = InitialConditionsFactory::template create<Mesh>
      (icType,paramsManager,reader,mpiCom);

    if (ic == static_cast<IC*>(NULL))
      {
	dice::glb::console->printFlush<dice::LOG_ERROR>
	  ("Invalid initial conditions type: '%s'.\n",icType.c_str());
	dice::glb::console->print<dice::LOG_ERROR>
	  ("init.type must be one of the following: %s.\n",
	   InitialConditionsFactory::getList().c_str());
	exit(-1);
      }
   
    typename IC::Params icp=ic->getParams();

    // Initialize mesh default params and parse them
    std::copy_n(&icp.x0.front(),NDIM,meshParams.x0);
    std::copy_n(&icp.delta.front(),NDIM,meshParams.delta);    
    meshParams.parse(reader,paramsManager);

    // setup cosmology if needed
    units.useCosmo = icp.useCosmo;
    units.useCosmo = paramsManager.
      get("cosmo",parserCategory(),units.useCosmo,reader,
	  PM::FILE_FIRST,"Use cosmological expansion ?");
    
    if (units.useCosmo)
      {	
	aStart = icp.aStart;
	aEnd = 1.0;
	//double dTau = 1.0;
	
	aStart = paramsManager.
	  get("aStart",parserCategory(),aStart,reader,
	      PM::FILE_FIRST,"Initial value of the scale factor");

	aEnd = paramsManager.
	  get("aEnd",parserCategory(),aEnd,reader,
	      PM::PARSER_FIRST,"Final value of the scale factor");
	
	double dLogA=0.1;
	dLogA = paramsManager.
	  get("dLogA",parserCategory(),dLogA,reader,
	      PM::PARSER_FIRST,"Maximum value for d(log(a))=da/a.");

	icp.cosmoParams.parse(reader,paramsManager);
	units.cosmology.initialize(icp.cosmoParams);
	
    	solverInterfaceParams.t0=units.cosmology.tau_of_a(aStart);
	solverInterfaceParams.tEnd=units.cosmology.tau_of_a(aEnd);
	solverInterfaceParams.dt=dLogA;
	
	solverInterfaceParams.parse(reader,paramsManager,true);

	units.H = (icp.H!=1)?(icp.H):1;
	units.H = paramsManager.
	  get("H","units",units.H,reader,
	      PM::FILE_FIRST,
	      "Hubble parameter (default value in s-1).");
      }
    else
      {	
	solverInterfaceParams.dt=1.0E-2; // set a default value 
	solverInterfaceParams.parse(reader,paramsManager);

	units.G = icp.G;
	units.G = paramsManager.
	  get("G","units",units.G,reader,
	      PM::FILE_FIRST,
	      "Gravitational constant (physical value: 6.67384E-11 m3 kg-1 s-2).");
      }

    char comment[256];
    sprintf(comment,"Length unit expressed in meters (default: %s).",
	    icp.defaultUnitLength.c_str());
    units.length = icp.unitLength;    
    units.length = paramsManager.
      get("length","units",units.length,reader,
	  PM::FILE_FIRST,comment);

    sprintf(comment,"Velocity unit expressed in meters (default: %s).",
	    icp.defaultUnitVelocity.c_str());
    units.velocity = icp.unitVelocity;
    units.velocity = paramsManager.
      get("velocity","units",units.velocity,reader,
	  PM::FILE_FIRST,comment);

    sprintf(comment,"Mass unit expressed in meters (default: %s).",
	    icp.defaultUnitMass.c_str());
    units.mass = icp.unitMass;
    units.mass = paramsManager.
      get("mass","units",units.mass,reader,
	  PM::FILE_FIRST,comment);
    
    units.update();

    mass = icp.mass;
    mass = paramsManager.
      get("mass",parserCategory(),mass,reader,
	  PM::FILE_FIRST,
	  "Initial total mass");

    fileDumps.setOutputDir(solverInterfaceParams.outputDir);

    dtMax = solverInterfaceParams.dt;

    symmetry = dice::RGSV::NONE;
    dice::RegularGridSymmetrySelect symmetryTypeSelect;   
    symmetryStr = paramsManager.
      get("symmetry",parserCategory(),
	  symmetryTypeSelect.getString(symmetry,true),reader,
	  PM::PARSER_FIRST,
	  symmetryTypeSelect.getAllString("Type of symmetry to enforce (%s)"));
    symmetry = symmetryTypeSelect.getVal(symmetryStr,true);

    cflCondition = CflConditionTypeV::CFL_RHOMAX;
    CflConditionTypeSelect cflConditionTypeSelect;
    std::string cflConditionStr = paramsManager.
      get("cflCondition",parserCategory(),
	  cflConditionTypeSelect.getString(cflCondition,true),reader,
	  PM::PARSER_FIRST,
	  cflConditionTypeSelect.getAllString("CFL condition type (%s)"));
    cflCondition = cflConditionTypeSelect.getVal(cflConditionStr,true);
    
    if ((cflCondition == CflConditionTypeV::RHOMAX)||
	(cflCondition == CflConditionTypeV::CFL_RHOMAX))
      {
	cflRhoMax = 0.01;//solverInterfaceParams.dt;
	cflRhoMax = paramsManager.
	  get("cflRhoMax",parserCategory(),cflRhoMax,reader,
	      PM::PARSER_FIRST,
	      "Value of the maximum density condition (fraction of an orbit size)");
      }
    
    if ((cflCondition == CflConditionTypeV::CFL)||
	(cflCondition == CflConditionTypeV::CFL_RHOMAX))
      {
	cflGrid = 0.25;
	cflGrid = paramsManager.
	  get("cflGrid",parserCategory(),cflGrid,reader,
	      PM::PARSER_FIRST,
	      "Value of the CFL condition (fraction of a FFT grid pixel size)");
      }
    
    fileDumps.parseFromManager(paramsManager,reader,classVersion(),serializedVersion);
    
    refineThreshold = 3.0;
    refineThreshold = paramsManager.
      get("refineThreshold",parserCategory(),refineThreshold,reader,
	  PM::PARSER_FIRST,
	  "NOT used");
    
    coarsenHysteresis = 0.95;
    coarsenHysteresis = paramsManager.
      get("coarsenHysteresis",parserCategory(),coarsenHysteresis,reader,
	  PM::PARSER_FIRST,
	  "NOT used"); 

    volumeThreshold = 1.E+10; // <=> desactivated
    volumeThreshold = paramsManager.
      get("volumeThreshold",parserCategory(),volumeThreshold,reader,
	  PM::PARSER_FIRST,
	  "Volume threshold to trigger refinement (NOT used)");

    invariantThreshold = 1.E-6;
    invariantThreshold = paramsManager.
      get("invariantThreshold",parserCategory(),invariantThreshold,reader,
	  PM::PARSER_FIRST,
	  "Poincare invariant threshold to trigger refinement");

    splitLongestEdge = SPLIT_LONGEST_EDGE_DEFAULT;
    splitLongestEdge = paramsManager.
      get("splitLongestEdge",parserCategory(),splitLongestEdge,reader,
	  PM::PARSER_FIRST,
	  "Set to split the longest longest edge when refining instead of trying to minimize the invariant");

    maxSimplexLevel = -1;
    maxSimplexLevel = paramsManager.
      get("maxSimplexLevel",parserCategory(),maxSimplexLevel,reader,
	  PM::PARSER_FIRST,
	  "Maximum level of refinement a simplex is allowed to reach (-1  for no maximum, unrefined simplices have level 0). See also 'projection.useVerticesThreshold'.",
	  serializedVersion>0.215);

    pUseVerticesThreshold = maxSimplexLevel;
    pUseVerticesThreshold = paramsManager.
      get("useVerticesThreshold","projection",pUseVerticesThreshold,reader,
	  PM::PARSER_FIRST,
	  "Minimum refinement level of a simplex above which its vertices are used to project it instead of exact integration. This should be equal to either '-1' or the value of 'solver.maxSimplexLevel'.",
	  serializedVersion>0.215);

    gatheredPotentialAllocLimit = 0;
    gatheredPotentialAllocLimit = paramsManager.
      get("gatheredPotentialAllocLimit",parserCategory(),gatheredPotentialAllocLimit,reader,
	  PM::PARSER_FIRST,
	  "Maximum (approximate) amount of memory in GigaBytes reserved for gathering the potential on local MPI processes. The larger, the faster, multiple passes being used to compensate for the lack of memory if needed. Set to 0 for unlimited, which amounts to allocating locally a grid equivalent to the FFT grid ( = 2^(NDIM*fftGridLevel)*sizeof(double) bytes). THIS OPTION IS DISABLED FOR NOW.",
	  serializedVersion>0.215);
    
    fftGridLevel = D_AMR_ROOT_LEVEL+2;
    fftGridLevel = paramsManager.
      get("fftGridLevel",parserCategory(),fftGridLevel,reader,
	  PM::PARSER_FIRST,
	  "Size of the FFT grid used for the Poisson solver (effective size is 2^fftGridLevel)");

    // Set the CFL condition max displacement SIZE
    cflSizeMax = cflGrid * meshParams.delta[0] / (1<<fftGridLevel);

    maxAmrLevel = fftGridLevel;
    maxAmrLevel = paramsManager.
      get("maxAmrLevel",parserCategory(),maxAmrLevel,reader,
	  PM::IGNORE_FILE,
	  "Maximum refinement level of the AMR grid (don't touch that, it's already set ;) )");
    
    rebuildAmrEvery = 1;    
    rebuildAmrEvery = paramsManager.
      get("rebuildAmrEvery",parserCategory(),rebuildAmrEvery,reader,
	  PM::PARSER_FIRST,
	  "How many timestep to wait before entirely rebuilding the AMR grids. (IGNORED : not implemented yet!)",
	  serializedVersion>0.105);
    
    double nThreads = dice::glb::num_omp_threads;
    nThreads = paramsManager.
      get("nThreads",parserCategory(),nThreads,reader,
	  PM::PARSER_FIRST,
	  "The number of threads to use per process. This parameter will override the default 'threads_per_node' parameter and will be saved from one run to another.",
	  serializedVersion>0.115);
    // The default number of threads for the entire library is overriden
    dice::setGlobalNumThreads(nThreads);
    
    fastAmrBuild = 1;    
    fastAmrBuild = paramsManager.
      get("fastAmrBuild",parserCategory(),fastAmrBuild,reader,
	  PM::PARSER_FIRST,
	  "If true, refine the AMR grid so that it only fits roughly the resolution of the mesh. This is probably always good enough and much faster ...",
	  serializedVersion>0.115);

    checkProjectedDensity = 1;
    checkProjectedDensity = paramsManager.
      get("checkProjectedDensity",parserCategory(),checkProjectedDensity,reader,
	  PM::PARSER_FIRST,
	  "Check the positivity of the projected density (also look for NaN values).",
	  serializedVersion>0.105);            

    phSortThreshold = 0.05;
    phSortThreshold = paramsManager.
      get("phSortThreshold",parserCategory(),phSortThreshold,reader,
	  PM::PARSER_FIRST,
	  "Maximum allowed ratio of randomly distributed to ordered cells before triggering a Peano-Hilbert sort of the local meshes.",
	  serializedVersion>0.115); 

    fftWisdom=1;
    fftWisdom=paramsManager.
      get("fftWisdom",parserCategory(),fftWisdom,reader,
	  PM::PARSER_FIRST,
	  "Which level of FFTWISDOM to use, in the range {0,1,2,3}={ESTIMATE,MEASURE,PATIENT,EXHAUSTIVE}. The higher the faster, but the slower the initialization ...",
	  serializedVersion>0.115);  
    
    projectionOrder=1;
    projectionOrder=paramsManager.
      get("projectionOrder",parserCategory(),projectionOrder,reader,
	  PM::PARSER_FIRST,
	  "Order of the exact mass projection (only 0 and 1 are implemented so far)",
	  serializedVersion>0.115);   

    exportFftWisdom=std::string("wisdom.fftw");
    exportFftWisdom=paramsManager.
      get("exportFftWisdom",parserCategory(),exportFftWisdom,reader,
	  PM::PARSER_FIRST,
	  "The name of a file to which the acquired FFTW wisdom will be exported. Set to 'SKIP' or 'skip' to prevent exportation.",
	  serializedVersion>0.135); 

    importFftWisdom=std::string("wisdom.fftw");
    importFftWisdom=paramsManager.
      get("importFftWisdom",parserCategory(),importFftWisdom,reader,
	  PM::PARSER_FIRST,
	  "The name of a file from which FFTW wisdom will be imported. Set to 'SKIP' or 'skip' to prevent importation.",
	  serializedVersion>0.135); 

    noRepartWeight=1;
    noRepartWeight=paramsManager.
      get("noRepartWeight",parserCategory(),noRepartWeight,reader,
	  PM::PARSER_FIRST,
	  "If true, MPI regions weight is directly proportional to the number of simplices.",
	  serializedVersion>0.175); 

    sprintf(comment,"Required accuracy level for projection %s",
	    D_ENABLE_ACCURACY_CHECKING?"(ENABLED).":"(DISABLED). Use D_ENABLE_ACCURACY_CHECKING compile time option to enable.");
    accuracyLevel = 1.E-3;
    accuracyLevel = paramsManager.
      get("accuracyLevel","projection",accuracyLevel,reader,
	  PM::PARSER_FIRST,
	  comment,
	  serializedVersion>0.185); 
    /*
    pForceHighResolution=-1;
    pForceHighResolution=paramsManager.
      get("forceHighResolution","projection",pForceHighResolution,reader,
	  PM::PARSER_FIRST,
	  "Forces the projection to be entirely computed using high resolution floating point type (" STRINGIFY(D_PROJECTION_HR_FLOAT_TYPE) ").",
	  serializedVersion>0.195);
    */    

    pEnableHRModeThreshold=0;
    pEnableHRModeThreshold=paramsManager.
      get("enableHRModeThreshold","projection",pEnableHRModeThreshold,reader,
	  PM::PARSER_FIRST,
	  "Fraction of the simplices that are allowed to be reprojected before definitively switching to high resolution projection mode (between 0 and 1, negative number=>always).",
	  serializedVersion>0.195);
    pForceHighResolution=(pEnableHRModeThreshold<0)?true:false;

    pResetHRMode=0;
    pEnableHRModeThreshold=paramsManager.
      get("resetHRMode","projection",pResetHRMode,reader,
	  PM::PARSER_FIRST,
	  "Ignore whether high resolution mode was set on a previous run and set it to false. This is only meaningfull when restarting a run.",
	  serializedVersion>0.195);

    /*
    pSwitchHRModeFileName="SWITCH_HRMODE";
    pSwitchHRModeFileName=paramsManager.
      get("switchHRModeFileName","projection",pSwitchHRModeFileName,reader,
	  PM::PARSER_FIRST,
	  "The name of the file to create in the output directory in order to toggle 'projection.forceHighResolution' on the fly.",
	  serializedVersion>0.195); 
    */
    pAnisotropyThreshold=1.e6;
    pAnisotropyThreshold=paramsManager.
      get("anisotropyThreshold","projection",pAnisotropyThreshold,reader,
	  PM::PARSER_FIRST,
	  "The simplex anisotropy threshold up to which exact projection is used (expressed as the maximum to minimum simplex extent ratio)",
	  serializedVersion>0.165); 

    pWidthThreshold=1.e-7;
    pWidthThreshold=paramsManager.
      get("widthThreshold","projection",pWidthThreshold,reader,
	  PM::PARSER_FIRST,
	  "The width of an element, expressed as a fraction of the bounding box size, down to which exact projection is used",
	  serializedVersion>0.165);    

    specifyProfileCenter=0;
    specifyProfileCenter=paramsManager.
      get("specifyProfileCenter",FileDumps::parserCategory(),specifyProfileCenter,reader,
	  PM::PARSER_FIRST,
	  "Profile is centered on the center of the bounding box if false, or use user specified coordinates if true (see densityProfileCenter)",
	  serializedVersion>0.205); 
    
    for (int i=0;i<NDIM;++i)
      {
	densityProfileCenter[i]=0;
	densityProfileCenter[i]=paramsManager.
	  get("densityProfileCenter",FileDumps::parserCategory(),densityProfileCenter[i],i,
	      reader,PM::PARSER_FIRST,
	      "The coodinates of the density profile center (only used if specifyProfileCenter is set)",
	      serializedVersion>0.205); 
      }

    // Non managed parameters (not saved in restart files)
    typename PM::Parser *parser = paramsManager.getParser();

    skipInitialPoisson=1;
    skipInitialPoisson=parser->
      get("skipInitialPoisson",parserCategory(),skipInitialPoisson,
	  "Skips the pre-solving of poisson equation for when the CFL conditon is set to CFL_RHOMAX. Instead, measure initial RHOMAX directly on the sheet.");
    
    squeeze=0;
    squeeze=parser->get("squeezeMesh",parserCategory(),squeeze,
			"You don't want to know ;)");

    dumpInitialMesh=0;
    dumpInitialMesh=parser->
      get("dumpInitialMesh",parserCategory(),dumpInitialMesh,
	  "Dump the mesh right after creating the initial conditions.");
  }
  
  bool checkCoarsen(std::vector<Simplex *> &s1, std::vector<Simplex *> &s2,
		    Vertex *v1, Vertex *v2, Vertex *v)
  {
    return geometry->template distance2<Coord,NDIM_W>
      (v1->getCoordsConstPtr(),v2->getCoordsConstPtr()) <= coarsenThreshold2;
  }
    
  double checkRefine_getValue(Simplex *s)
  {
    double result=0;

    if ((maxSimplexLevel<0)||(s->getLevel()<maxSimplexLevel))
      {
	result = dice::slv::refine::
	  poincareInvariantWithSegTracers_order1<Mesh>(s,geometry).first;
	
	result *= invariantThreshold_inv;
	if (result < 1) result=0;    
      }
    
    return result;
  }
    
  //CheckRefineReturnType 
  int checkRefine_getSplitSegmentIndex(Simplex *s, double &oldVal)
  {   
    int ret;

    // invariant is above the threshold, this simplex should be refined !
    if (splitLongestEdge)
      {
	// Just refine the longest edge ...
	//ret = dice::slv::refine::length2<Mesh>(s,geometry).second; // Eulerian 
	ret=dice::slv::refine::lagrangianLength2<Mesh>(s,geometry).second; // Lagrangian
      }
    else
      {
	// The function we need to minimize when splitting
	struct Functor
	{
	  Functor(Mesh *m)
	  {
	    geometry=m->getGeometry();
	  }

	  double operator()(Simplex *s) const
	  {
	    return dice::slv::refine::
	      poincareInvariantWithSegTracers_order1<Mesh>
	      (s,geometry).first;
	  }
	      
	  GeometricProperties *geometry;
	} functor(mesh);
	    
	// Select the edge that minimize refinement criterion
	ret = dice::slv::refine::
	  findSplitSegment(mesh,s,functor,invariantThreshold);
	    
	// This will happen if splitting the simplex does not improve the invariant
	if (ret<0)
	  {	
	    ret=dice::slv::refine::lagrangianLength2<Mesh>(s,geometry).second;
	    /*
	      dice::glb::console->print<dice::LOG_STD_ALL>
	      ("Splitting longest edge @%ld (linear: %g, quadratic : %g)\n",
	      (long)s->getLocalIndex(),oldVal,
	      dice::slv::refine::poincareInvariantWithSegTracers<Mesh>(s,geometry)
	      );	
	    */
	  }
	   
      }
   
    return ret;
  }

  template <class W>
  void onWrite(W *writer)
  {
    writer->write(&densityMax);
    writer->write(&initialMeshResolution,NDIM);
    writer->write(&potentialEnergy);
    writer->write(&expansionEnergy);
    writer->write(&curSolverDeltaT);
    writer->write(&pForceHighResolution);
    fileDumps.write(writer);    
  }

  template <class R>
  void onRead(R *reader)
  {    
    reader->read(&densityMax);
    reader->read(&initialMeshResolution,NDIM);
    reader->read(&potentialEnergy);
    reader->read(&expansionEnergy);  
    reader->read(&curSolverDeltaT);

    int savedValue=pForceHighResolution;
    if (serializedVersion>0.195) 
      reader->read(&pForceHighResolution);
    if (pResetHRMode) pForceHighResolution=savedValue;

    fileDumps.read(reader);    
  }

  void onBuildMesh(Mesh *m, const MeshParams &meshParams)
  {
    auto *implicitTesselation=ic->createImplicitTesselation();
    m->build(implicitTesselation,meshParams);
    implicitTesselation->getResolution(initialMeshResolution);
    delete implicitTesselation;
  }
   
  template <class SP>
  void onInitialize(Mesh *m,const SP& solverParams, double t, bool restart)
  {          
    init_impl(m,solverParams,t,restart);
    sortMesh();
  }

  template <class OutputStream>
  void onWriteTimings(double t, OutputStream &out, bool writeHeader)
  {
    if (writeHeader)
      {
	if (units.useCosmo) out << " " << "a";
      }
    else
      {
	if (units.useCosmo) out << " " << units.cosmology.a_of_tau(t,aStart);
      }
  }

  template <class OutputStream>
  void onWriteStatistics(double t, OutputStream &out, bool writeHeader)
  {
 
    if (writeHeader)
      {
	if (units.useCosmo) out << " " << "a";
	out << " " << "pMassErr"
	    << " " << "V"
	    << " " << "V_Q"
	    << " " << "Ek"
	    << " " << "Ek_Q"
	    << " " << "Ep"
	    << " " << "E"
	    << " " << "E_Q";
      }
    else
      {
	if (units.useCosmo) out << " " << units.cosmology.a_of_tau(t,aStart);
	out << " " << projectedMassErr
	    << " " << totalVolume
	    << " " << totalVolume_Q
	    << " " << kineticEnergy
	    << " " << kineticEnergy_Q
	    << " " << potentialEnergy
	    << " " << potentialEnergy+kineticEnergy+expansionEnergy
	    << " " << potentialEnergy+kineticEnergy_Q+expansionEnergy;
      }    
  }

  void onResimulate()
  {    
    FOREACH_THREAD_VERTEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
      {	
	for(;it!=it_end;++it)
	  {
	    Vertex *v=(*it);
	    Coord *c=v->getCoordsPtr();
	    Coord *ic=v->initCoords.getPointer();
		
	    for (int k=0;k<NDIM_W;++k) c[k]=ic[k];		  
	  }
      }
    FOREACH_THREAD_END;
    
    initCellData(true);
  }

  void onNewTimeStep(long stepId, double curTime, double &newDeltaT)
  { 
    updateTimeStep(curTime,newDeltaT);
    updateInvariantThreshold(curTime);    

    curStepIndex=stepId;
    curSolverTime=curTime;

    oldSolverDeltaT=curSolverDeltaT;
    curSolverDeltaT=newDeltaT;
       
    fileDumps.updateState(curStepIndex,curSolverTime,curSolverDeltaT,mpiCom);
  }

  void onAdvance()
  {       
    advance_impl();
  }

  void afterCoarsen(long nCoarsened) 
  {}
 
  void afterRefine(long nRefined) 
  {  
    newMeshSimplicesCount += nRefined;
  }
  
  // This function should return the TOTAL weight for the local region
  // For default weighting, return 0.
  // Set force to true to enforce repartitionning
  double onCheckRepartLocalWeight(double stepDuration, bool &force)
  {
    if (noRepartWeight) return 0;

    double waitDuration=projectBarrierTimer->lastSpent();
    double projectDuration=projectTimer->lastSpent();
    double otherDuration=stepDuration-waitDuration-projectDuration;    

    double otherWeight=mpiCom->min(otherDuration/mesh->getNSimplices());    
    double projectWeight=projectDuration/mesh->getNSimplices();
    
    // We may have (otherWeight<0) when using e.g. static potential solver
    double result = (otherWeight<0)?0:(otherWeight + projectWeight)*mesh->getNSimplices();
    
    return result;

    //return stepDuration - projectBarrierWaitDuration;
    //return -1;
  }

  // status is true if the GLOBAL mesh was repartitionned 
  void afterRepart(bool status)
  {
    repartStatus=status;
    
    // If the mesh was repartitionned, we should always sort it locally ...
    if (status) sortMesh();
    else 
      {
	// We also sort the mesh if the fraction of unsorted cells is higher
	// than the threshold
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Checking local mesh ordering ... ");
	double ratio = mpiCom->max(double(newMeshSimplicesCount)/mesh->getNSimplices());
	if (ratio>phSortThreshold)
	  {
	    dice::glb::console->print<dice::LOG_STD>
	      ("threshold crossed (f=%.2g).\n",ratio);
	    sortMesh();
	  
	    newMeshSimplicesCount=0;
	    //ratio = mpiCom->max(double(newMeshSimplicesCount)/mesh->getNSimplices());
	  }
	else
	  dice::glb::console->print<dice::LOG_STD>
	    ("done. (f=%.2g < %.2g)\n",ratio,phSortThreshold);	
      }

  }
  
 
protected:

  template <class SP>
  void init_impl(Mesh *m, const SP& solverParams, double t, bool restart,
		 int forceNoFFTLevel=0)
  {
    mesh = m;
    geometry = mesh->getGeometry();        

    const MeshParams &p=mesh->getParams();

    // AMR grid initialization
    std::vector<double> x0;
    std::vector<double> delta;

    mesh->getBoundingBox(std::back_inserter(x0),std::back_inserter(delta)); 
 
    localAmrDensity.initialize(x0.begin(),delta.begin());
    // Potential grid initialization
    // setup grid parameters
    RegularGridParams rgp;   
    std::copy(x0.begin(),x0.end(),rgp.x0);
    std::copy(delta.begin(),delta.end(),rgp.delta); 
    std::fill(rgp.resolution,rgp.resolution+NDIM,(1<<fftGridLevel));    

    if (forceNoFFTLevel==0) 
      {
	if ((importFftWisdom != std::string("SKIP"))&&
	(importFftWisdom != std::string("skip")))
	  fftSolver.importWisdom(importFftWisdom.c_str());
	else
	  dice::glb::console->print<dice::LOG_STD>("Importing FFTW wisdom: SKIPPED.\n");

	// Let the FFT solver slice the density grid for us ...
	const double pi=4.0*atan(1.0);
	FFTPoissonKernel kernel((units.useCosmo)?-1.0:-4.0*pi*units.G);
	fftSolver.initializeGridAndPlans(&potential,rgp,kernel,fftWisdom,false);  
	potential.clone(clonedDensity,false);

	if ((exportFftWisdom != std::string("SKIP"))&&
	    (exportFftWisdom != std::string("skip")))
	  fftSolver.exportWisdom(exportFftWisdom.c_str());
	else
	  dice::glb::console->print<dice::LOG_STD>("Exporting FFTW wisdom: SKIPPED.\n");
      }
    else if (forceNoFFTLevel>0)
      {
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Generating density grid using default slicer ... ");
	std::fill(rgp.resolution,rgp.resolution+NDIM,(1<<forceNoFFTLevel));
	typename RegularGrid::DefaultSlicer slicer(rgp,mpiCom->size(),mpiCom,
						   dice::glb::num_omp_threads,
						   Periodic);
	potential.initialize(slicer);
	
	dice::glb::console->printFlush<dice::LOG_STD>("done.\n");
      }
    
    // Setup some global parameters
    updateInvariantThreshold(t);
    refineThreshold2 = refineThreshold * (p.delta[0]/(1<<fftGridLevel) * sqrt(NDIM));
    refineThreshold2 = pow(refineThreshold2,2.0); 
    invariantThreshold_inv = 1.0F / (invariantThreshold);
    volumeThreshold_inv = 1.0F / (volumeThreshold);
    maximumSegmentLength2_inv = 1.0F / (refineThreshold2);
    coarsenThreshold2 = refineThreshold2*pow(coarsenHysteresis,2.0);
    newMeshSimplicesCount=0;

    // Initialize timers
    amrBuildTimer  = dice::glb::timerPool->pop("amrBuild");
    updateDensityTimer  = dice::glb::timerPool->pop("updateDensity");    
    projectTimer   = dice::glb::timerPool->pop("project");
    projectBarrierTimer = dice::glb::timerPool->pop("project_barrier");
    kickAndDriftTimer = dice::glb::timerPool->pop("kickDrift");
    driftTimer = dice::glb::timerPool->pop("drift");    
    fftSolverTimer = dice::glb::timerPool->pop("FFT");
    solvePoissonTimer = dice::glb::timerPool->pop("poisson");
    scatterDensityTimer = dice::glb::timerPool->pop("scatter_density");
    gatherPotentialTimer = dice::glb::timerPool->pop("gather_potential");
    sortMeshTimer = dice::glb::timerPool->pop("sort_mesh");
    statisticsTimer = dice::glb::timerPool->pop("stats");
    dumpTimer = dice::glb::timerPool->pop("dump");

    //mesh->dumpToNDnetwork((std::string("INIT_")+dumpManager.meshName).c_str());

    // Add some data functors to the mesh (quantites that can be computed on the 
    // fly, so they can be requested at any time but do not use extra storage)
    mesh->template 
      addCellDataFunctor< dice::slv::refine::
			  CellDataFunctor_poincareInvariantWithSegTracers_order1T >();
    mesh->template 
      addCellDataFunctor< ProjectedDensityFunctorT >();
    
    /*
      // These may be useful for debugging ...
    mesh->template 
      addCellDataFunctor< dice::slv::refine::
			  CellDataFunctor_poincareInvariantT >();
    mesh->template 
      addCellDataFunctor< dice::slv::refine::
			  CellDataFunctor_poincareInvariantWithSegTracersT >();
    mesh->template 
      addCellDataFunctor< dice::slv::refine::CellDataFunctor_lengthT >();
    mesh->template 
      addCellDataFunctor< dice::slv::refine::
			  CellDataFunctor_poincareInvariantBelowSegTracersCornerT >();
    mesh->template 
      addCellDataFunctor< dice::slv::refine::CellDataFunctor_volumeBelowT >();   
    mesh->template 
      addCellDataFunctor< SmoothedProjectedDensityFunctorT >();   
    */
    
    // Do that only if we are not restarting from a run
    if (!restart) 
      {	
	expansionEnergy=0;
	curSolverDeltaT=0;
	oldSolverDeltaT=0;
	potentialEnergy=0;
	
	initCellData(t);
	if (squeeze) squeezeMesh();  

	if ((cflCondition == CflConditionTypeV::CFL_RHOMAX)||
	    (cflCondition == CflConditionTypeV::RHOMAX))
	  {
	    if (!skipInitialPoisson) solvePoisson(t);
	    else
	      {
		std::vector<double> result(dice::glb::num_omp_threads,0);

		FOREACH_THREAD_SIMPLEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
		  {
		    const SimplexFunctor &projectedDensityFunctor=
		      *mesh->getSimplexFunctorPtr("projectedDensity");
		    
		    for (;it!=it_end;++it)
		    {
		      double val=projectedDensityFunctor.get(*it);
		      if (val>result[th]) result[th]=val;
		    }
		  }
		FOREACH_THREAD_END
	
		densityMax=result[0];
		for (int th=1;th<dice::glb::num_omp_threads;th++)
		  if (result[th]>densityMax) densityMax=result[th];

		if (units.useCosmo)
		  {
		    double a=units.cosmology.a_of_tau(t,aStart);
		    double fact = (3.0*a*units.cosmology.getParams().oM)/2.0;
		    double volume=1.0L;
		    for (int i=0;i<NDIM;++i) volume *= mesh->getParams().delta[i];
		    double avg = mass/volume;
		    densityMax = (densityMax/avg) * fact;
		  }
	      }
	  }
      }
    
    if (dumpInitialMesh)
      mesh->dumpToNDnetwork((fileDumps.getLocalFName(FileDumps::Mesh)+
			     std::string("_init")).c_str());    
  }

  void initCellData(double t, bool resimulate=false)
  {    
     std::vector<double> x0;
    std::vector<double> delta;

    mesh->getBoundingBox(std::back_inserter(x0),std::back_inserter(delta));     
    
    const long nPasses=ic->getNPasses();

    dice::glb::console->printFlush<dice::LOG_STD>
      ("Initializing cells (%s) in %ld pass(es) ... ",ic->getName().c_str(),nPasses);
    
    for (long pass=0; pass<nPasses; ++pass)
      {

	ic->initialize(&x0[0],&delta[0],t,units,pass);
	
	if (pass>0)
	  dice::glb::console->printFlush<dice::LOG_STD>("(%d%%)",(100*pass)/nPasses);    
	  	
	// First initialize simplex data
	FOREACH_THREAD_SIMPLEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
	  {	
	    for (;it!=it_end;++it)
	      {	    
		Simplex *s=(*it);
	    
		// Initialize mass before displacing the particles
		if (pass==0) 
		  {
		    double initialDensity=1;
		    // Mass is initially proportional to the lagrangian volume !
		    s->mass.init(mesh,s,&initialDensity);
		  }

		// Displace the tracer !
#ifndef NO_SIMPLEX_TRACERS
		if (pass==0) s->tracer.init(mesh,s); 
		typename Simplex::Tracer::Type *c = s->tracer.getPointer();	   
		ic->displace(c,s);
		
		geometry->sanitizeBoundary(c);
#endif
	    
		// and the seg tracers
		if (pass==0) s->segTracers.init(mesh,s);
		typename Simplex::SegTracers::Type *cst = s->segTracers.getPointer();
		for (int j=0;j<Simplex::NSEG;++j)
		  {
		    typename Simplex::SegTracers::Type *c = &cst[j*NDIM_W];	
		    ic->displace(c,s);
		    
		    geometry->sanitizeBoundary(c);
		  }
	      }
	  }
	FOREACH_THREAD_END;

	// Then vertices data and coordinates
	// note that projectedDensity temporarily holds the actual simplex density
	FOREACH_THREAD_VERTEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
	  {		
	    for (;it!=it_end;++it)
	      {
		Vertex *v=(*it);
		Coord *c=v->getCoordsPtr();

		if ((!resimulate)&&(pass==0)) 
		  {
		    v->initCoords.init(mesh,v);		    
		    const double defaultDensity=1;
		    v->projectedDensity.init(mesh,v,&defaultDensity);
		  }

		if (!resimulate)
		  ic->displace(c,*v->projectedDensity.getPointer(),v);
		else
		  std::copy_n(v->initCoords.getPointer(),NDIM_W,c);
		
		geometry->sanitizeBoundary(c);
	      }	
	  } 
	FOREACH_THREAD_END;
      }    
    ic->release();

    // Now if the mass is determined from the density set at initial conditions, we need 
    // to initialize it for each simplex from the vertices densities
    if (!resimulate)
      {
	double massRatio=0;

	FOREACH_THREAD_SIMPLEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
	  {
	    double m=0;
	    for (;it!=it_end;++it)
	      {
		Simplex *s=(*it);
		
		double density=s->getVertex(0)->projectedDensity.getValue();
		for (int i=1;i<Simplex::NVERT;++i)
		  density+=s->getVertex(i)->projectedDensity.getValue();
		density/=Simplex::NVERT;

		if (ic->useLagrangianMass())
		  {
		    s->mass.setValue(s->mass.getValue()*density);
		  }
		else
		  {
		    s->mass.init(mesh,s,&density);
		  }
		
		if (s->isLocal()) m+=s->mass.getValue();
	      }
	
#pragma omp atomic update
	    massRatio+=m;
	  }
	FOREACH_THREAD_END;

	// normalize to the total mass
	massRatio = mass/mpiCom->sum(massRatio);
	FOREACH_THREAD_SIMPLEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
	  {
	    for (;it!=it_end;++it)
	      {
		Simplex *s=(*it);
		s->mass.setValue(s->mass.getValue()*massRatio);	   
	      }
	  }
	FOREACH_THREAD_END;
      }

    // Finally, compute projected density from simplices mass
    updateProjectedDensityField();

    if (dice::glb::debug)
      mesh->dumpToNDnetwork((fileDumps.getLocalFName(FileDumps::Mesh)+
			     std::string("_init")).c_str());
    
    dice::glb::console->print<dice::LOG_STD>(" done.\n");
  }
  
  // Sort the local mesh along a peano hilbert curve
  void sortMesh()
  {       
    typedef typename dice::PeanoHilbertT<NDIM> PH;
    double x0[NDIM];
    double deltaInv[NDIM];
    mesh->getBoundingBox(x0,deltaInv,false); 
    
    for (int i=0;i<NDIM;++i) 
      deltaInv[i]=1.0/deltaInv[i];
    
    // A functor that takes a simplex and returns it lagrangian coordinates along 
    // a peano hilbert curve
    auto sf=[x0,deltaInv](Simplex *s) -> typename PH::HCode
      {
	typename PH::HCode result;

	// Coords are the barycenter of lagrangian vertices coordinates
	Coord c[NDIM]={0};
	for (int i=0;i<Simplex::NVERT;++i)
	  {
	    const Coord *vc=s->getVertex(i)->initCoords.getPointer();
	    for (int j=0;j<NDIM;++j)
	      c[j] += vc[j];
	  }
	    	    
	for (int i=0;i<NDIM;++i)
	  c[i]=c[i]/Simplex::NVERT;
		
	PH::coordsToLength(c,result,x0,deltaInv);
	return result;
      };

    // A functor that takes a vertex and returns it coordinates along a peano 
    // hilbert curve
    auto vf=[x0,deltaInv](Vertex *v) -> typename PH::HCode
      {
	typename PH::HCode result;
	PH::coordsToLength(v->initCoords.getPointer(),result,x0,deltaInv);
	return result;	
      };

    dice::glb::console->print<dice::LOG_STD>
      ("Sorting mesh cells (Peano-Hilbert):\n");
    dice::glb::console->indent();
    sortMeshTimer->start();

    mesh->sort(sf,vf);

    double t=sortMeshTimer->stop();
    dice::glb::console->unIndent();
    dice::glb::console->print<dice::LOG_STD>
      ("Done in %.2gs.\n",t);
    
  }

  void squeezeMesh(double factor = 0.999L)
  {
    FOREACH_THREAD_VERTEX_LGS(mesh,dice::glb::num_omp_threads,th,it)
      {
	for (;it!=it_end;++it)
	  {
	    Coord *c=it->getCoordsPtr();
	    for (int i=0;i<NDIM;++i)
	      c[i]*=factor;	   
	  }
      }
    FOREACH_THREAD_END;
  }
 
protected:

  void buildAmrGrid(int amrLevel)
  {
    dice::glb::console->printFlush<dice::LOG_STD>
      ("Building adaptive grid from mesh (LR:%d/%d) ... ",(int)D_AMR_ROOT_LEVEL,amrLevel);
    amrBuildTimer->start();
    static int nStepsSinceLastRebuild=0;
    int level=amrLevel-LocalAmrGrid::ROOT_LEVEL;
    rebuildAmrEvery = 1; // ignore rebuildAmrEvery, not implemented ....

    double epsilon=localAmrDensity.getBBoxSize(0);
    for (int i=1;i<NDIM;++i)
      if (localAmrDensity.getBBoxSize(i)>epsilon)
	epsilon=localAmrDensity.getBBoxSize(i);
    epsilon*=pWidthThreshold;  
    
    double volumeThreshold=pow(epsilon,double(NDIM));
    double lengthThreshold=epsilon;
    
    dice::ProjectionTag projectionTag
      (pUseVerticesThreshold,dice::ProjectionTag::sampleFromVertices,
       volumeThreshold,dice::ProjectionTag::sampleWithOverlap,
       lengthThreshold,dice::ProjectionTag::sampleWithOverlap,
       pAnisotropyThreshold,dice::ProjectionTag::sampleWithOverlap);

    if ((repartStatus)||((rebuildAmrEvery-1)<=nStepsSinceLastRebuild))
      {
	dice::glb::console->printFlush<dice::LOG_STD>("(from scratch) ");
	if (fastAmrBuild)
	  {
	    dice::glb::console->printFlush<dice::LOG_STD>("(fast) ");
	    localAmrDensity.buildFromMesh_Fast
	      (mesh,projectionTag,1.0,level,true,dice::glb::num_omp_threads,true,NDIM>2);
	  }
	else
	  {
	    localAmrDensity.buildFromMesh
	      (mesh,1.0,level,true,dice::glb::num_omp_threads);
	  }
	nStepsSinceLastRebuild=0;
      }
    else 
      {
	// THIS IS NOT WORKING YET :
	// AMR grid update is not implemented  ...
	dice::PRINT_SRC_INFO(dice::LOG_ERROR);
	dice::glb::console->print<dice::LOG_ERROR>
	  ("rebuildAmrEvery must be 1 (NOT IMPLEMENTED YET).\n");
	exit(-1);
	/*
	//localAmrDensity.buildFromMesh(mesh,1.0,level,false);
	dice::localAmrGridVisitors::SetValueT<LocalAmrGrid> resetVisitor(0);
	localAmrDensity.visitTree(resetVisitor);
	*/
	nStepsSinceLastRebuild++;
      }
    double elapsed=amrBuildTimer->stop();
    dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed);     
  }

  template <typename IT =IFloat, 
	    typename HT =HFloat, 
	    typename ST =SFloat,
	    bool checkAccuracy=D_ENABLE_ACCURACY_CHECKING,
	    class LOG   =dice::LOG_WARNING
	    >
  int projectMesh(int pass=0)
  {
    long reprojectedSimplicesCount=0;

    projectTimer->start(); 
    if (projectionOrder==0)
      {
	const SimplexFunctor &projectedDensityFunctor=
	  *mesh->getSimplexFunctorPtr("projectedDensity");  

	reprojectedSimplicesCount=localAmrDensity.
	  template projectMesh0<Mesh,SimplexFunctor,
				checkAccuracy,IT,HT,ST>
	  (mesh,projectedDensityFunctor,accuracyLevel,
	   (fastAmrBuild)?true:false, // Use simplex tags
	   dice::glb::num_omp_threads,
	   true); 	   
      }
    else
      {
	const VertexFunctor &projectedDensityFunctor=
	  *mesh->getVertexFunctorPtr("projectedDensity");
	const SimplexFunctor &projectedDensityGradientFunctor=
	  *mesh->getSimplexFunctorPtr("projectedDensityGradient");  
	
	reprojectedSimplicesCount=localAmrDensity.
	  template projectMesh1<Mesh,VertexFunctor,SimplexFunctor,
				checkAccuracy,IT,HT,ST>
	  (mesh,projectedDensityFunctor,projectedDensityGradientFunctor,
	   accuracyLevel,
	   (fastAmrBuild)?true:false, // Use simplex tags
	   dice::glb::num_omp_threads,
	   true); 
      }
    double elapsed=projectTimer->stop();

    // Synchronize MPI processes if necessary
    // This gives us a measure of the projection imbalance.
    if (mpiCom->size() > 1)
      {
	projectBarrierTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>("Synchronizing MPI processes ... ");
	
	// This will force a synchronisation anyway, so no need for the barrier
	auto mms = mpiCom->minMaxSum<double,true>(elapsed);
	mpiCom->sum(reprojectedSimplicesCount);

	double min = mms.first.first;
	double max = mms.first.second;
	double avg = mms.second / mpiCom->size();
	double imbalance = max / avg;

	double projectBarrierWaitDuration = projectBarrierTimer->stop();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("done in %.3lgs.\n",projectBarrierWaitDuration);
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Global projection imbalance factor: %.3lg (Tmin=%.3lgs, Tmax=%.3lgs, Tavg=%.3lgs).\n",imbalance,min,max,avg);	
      }

    if (reprojectedSimplicesCount > pEnableHRModeThreshold * mesh->getGlobalNCells(NDIM))
      {
	dice::glb::console->print<dice::LOG_WARNING>
	  ("Reprojected simplices threshold reached, switching high resolution projection mode ON.\n");	
	pForceHighResolution=!pForceHighResolution;
      }

    if (checkProjectedDensity)
      {
	double threshold = -1.E-5;
	typename dice::TimerPool::Timer timer;
	timer.start();
	dice::glb::console->printFlush<dice::LOG_STD>("Checking projected mass ... ");
	typedef LessOrNanT<typename LocalAmrGrid::Data> LessOrNan;

	LessOrNan lessOrNan(threshold);
	dice::localAmrGridVisitors::CountIfT<LocalAmrGrid,LessOrNan> 
	  negativeVoxelsVisitor(lessOrNan);
	localAmrDensity.visitTree(negativeVoxelsVisitor);
	long count = negativeVoxelsVisitor.getCount();
	double elapsed=timer.stop();
	if (mpiCom->max(count)>0)
	  //if (count>0)
	  {	    	    
	    dice::glb::console->printFlush<dice::LOG_STD>("failed!\n");
	    if (count>0)
	      dice::glb::console->print<LOG>
		("Negative or NaN values found for the projected mass in %ld voxel(s).\n",
		 negativeVoxelsVisitor.getCount());
	    
	    dice::localAmrGridVisitors::GetVoxelIfT<LocalAmrGrid,LessOrNan> 
	      negativeGetVoxelsVisitor(lessOrNan);
	    localAmrDensity.visitTree(negativeGetVoxelsVisitor);
	    typedef typename LocalAmrGrid::Voxel Voxel;
	    const std::vector<Voxel*> &result = negativeGetVoxelsVisitor.getResult();
	    
	    if (result.size()>0)
	      {
		Coord coords[3]={0};
		localAmrDensity.index2CenterCoords(result[0]->getIndex(),&coords[0]);
		dice::glb::console->printFlush<LOG>
		  ("Negative voxels ID: V=%lg@%ld(%lg,%lg,%lg)",
		   result[0]->data,
		   result[0]->getIndex(),
		   coords[0],coords[1],coords[2]);
		for (long i=1;i<result.size();++i)
		  {
		    localAmrDensity.index2CenterCoords(result[i]->getIndex(),&coords[0]);
		    dice::glb::console->printFlush<LOG>
		      (", V=%lg@%ld(%lg,%lg,%lg)",
		       result[i]->data,
		       result[i]->getIndex(),
		       coords[0],coords[1],coords[2]);
		  }
		dice::glb::console->printFlush<LOG>(".\n");
	      }

	    if (pass==0)
	      {
		if (count>0)
		  dice::glb::console->print<LOG>
		    ("This most likely means that the floating point number type you are using for projection is not precise enough.\n");

		if (!dice::hlp::SameType<IT,HT>::value)
		  {
		    if (count>0)
		      dice::glb::console->print<LOG>
			("I Will try again using higher precision.\n");
			
		    mpiCom->barrier();

		    dice::glb::console->printFlush<dice::LOG_STD>
		      ("Retrying to project with high precision (expect up to a ~10 folds slowdown).\n");
		    // Recompute density with higher resolution
		    dice::glb::console->printFlush<dice::LOG_STD>
		      ("Updating projected density (HR) ... ");
		    updateProjectedDensityField<HT>();
		    dice::glb::console->printFlush<dice::LOG_STD>("done.\n");
		    // We need to rebuild the grid as the simplices cache has been erased
		    buildAmrGrid(maxAmrLevel);	
		
		    // and reproject using quad precision floats
		    return projectMesh<HT,HT,HT,false,dice::LOG_ERROR>(pass+1);
		  }
	      }
	    else 
	      {
		if (count>0)
		  dice::glb::console->print<LOG>
		    ("This may either mean that you need extremely high precision, or that there is a subtle bug in the projection. Better luck next time ;(\n");
		
		if (dice::glb::debug)
		  {
		    localAmrDensity.toVtk((fileDumps.getLocalFName(FileDumps::Amr)+
					   std::string("_dbg")).c_str());
		    mesh->dumpToNDnetwork((fileDumps.getLocalFName(FileDumps::Mesh)+
					   std::string("_dbg")).c_str());
		  }
		mpiCom->barrier();
		//exit(-1);
	      }	    	   	      
	  }
	dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);
      }    
    return pass;
  }
 
  void solvePoisson(double t)
  {    
    double elapsed;    
    static int counter=0;
    double a = t;
    double unitsK=units.mass * pow(units.velocity,2.0);
    double unitsP=units.mass * pow(units.length,2.0);

    if (units.useCosmo)
      {
	a=units.cosmology.a_of_tau(t,aStart);	
	unitsK*=pow(a/units.H,2.0)*units.cosmology.getParams().oM;
	unitsP*=units.cosmology.getParams().oM;//pow(units.length,NDIM);
	dice::glb::console->print<dice::LOG_STD>("Solving poisson equation (a=%g):\n",a);
      }
    else dice::glb::console->print<dice::LOG_STD>("Solving poisson equation:\n");

    dice::glb::console->indent();
    solvePoissonTimer->start();

    // Compute projected density and its gradient over the mesh
    dice::glb::console->printFlush<dice::LOG_STD>("Updating projected density ... ");
    updateDensityTimer->start();
    updateProjectedDensityField();
    elapsed=updateDensityTimer->stop();
    dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed); 
    
    // Compute statistics 
    dice::glb::console->printFlush<dice::LOG_STD>
      ("Computing mesh statistics ... ");
    statisticsTimer->start();
    
    double meshQuadratureResult[4];
    auto quadratures = dice::hlp::makeObjectList(KineticEnergyMQT<Mesh>(), 
						 KineticEnergyMQT<Mesh,2>(),
						 VolumeMQT<Mesh>(),
						 VolumeMQT<Mesh,2>());

    mesh->quadratureFromList
      (quadratures,meshQuadratureResult, dice::glb::num_omp_threads, true);
   
    
    kineticEnergy=meshQuadratureResult[0] * unitsK;
    kineticEnergy_Q=meshQuadratureResult[1] * unitsK;
    totalVolume=meshQuadratureResult[2];
    totalVolume_Q=meshQuadratureResult[3];

    elapsed = statisticsTimer->stop();
    dice::glb::console->print<dice::LOG_STD>("done in %lgs\n",elapsed);
    
    dice::glb::console->indent();
    dice::glb::console->print<dice::LOG_INFO>
      ("-> Kinetic energy Ek = %e / %e\n",kineticEnergy,kineticEnergy_Q);
    dice::glb::console->print<dice::LOG_INFO>
      ("-> Total volume Vt = %e / %e\n",totalVolume,totalVolume_Q);
    dice::glb::console->unIndent();

    dumpMeshDensityProfile();

    // For debugging purpose
    /*
    dice::PRINT_SRC_INFO(dice::LOG_WARNING);   
    printf("COMMENT ME HERE !!!!\n");
#pragma omp parallel for
    for (int j=0;j<dice::glb::num_omp_threads;j++)
      {	
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-10.0,10.0);

	const vertexPtr_LGS_iterator itv_end=mesh->vertexLGSEnd();
	for (vertexPtr_LGS_iterator it=
	       mesh->vertexLGSBegin(j,dice::glb::num_omp_threads);it!=itv_end;++it)
	  {
	    Vertex *v=(*it);
	    Coord *c=v->getCoordsPtr();
	    
	    if (v->getLocalIndex()==1212288)
	      {
		printf("%d: %.20e %.20e %.20e\n",v->getLocalIndex(),c[0],c[1],c[2]);
		c[0]=1.e-15;
		printf("%d: %.20e %.20e %.20e\n",v->getLocalIndex(),c[0],c[1],c[2]);
	      }
	  }
      }
    */
    // Refine the AMR grid to the resolution of the local projected mesh
    buildAmrGrid(maxAmrLevel);   
 
    /*
    // Force switching to high resolution projection via file signal
    // This requires MPI communications and is not usefull in general, kept
    // for debugging purpose
    if (mpiCom->rank()==0)
      {
	std::string fName=fileDumps.getOutputDir()+pSwitchHRModeFileName;
	if (dice::myIO::exists(fName)) 
	  {
	    dice::glb::console->print<dice::LOG_WARNING>
	      ("High resolution projection mode switch signal received.\n");	
	    dice::glb::console->print<dice::LOG_WARNING>
	      ("Switching high resolution projection mode %s.\n",
	       pForceHighResolution?"OFF":"ON");
	    pForceHighResolution=!pForceHighResolution;
	    
	    dice::myIO::rmFile(fName);
	  }
      }
    // Set a coherent projection mode on all mpi processes
    mpiCom->Bcast(&pForceHighResolution,0,1);
    */

    // Project the local mesh's 'projectedDensity' field onto the AMR grid
    if (pForceHighResolution)
      projectMesh<HFloat,HFloat,HFloat,false>();
    else
      projectMesh<IFloat,HFloat,SFloat>();
   
    if (fileDumps.checkEvent(FileDumps::Amr,true))
      {
	dumpTimer->start();
	localAmrDensity.toVtk(fileDumps.getLocalFName(FileDumps::Amr).c_str());
	double elapsed = dumpTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>
	  ("File was dumped in %lgs.\n",elapsed);
      }

    // Convert the AMR local grids to a single MPI distributed regular grid
    scatterDensityTimer->start();
    dice::glb::console->printFlush<dice::LOG_STD>
      ("Scattering density field (AMR->grid):\n");


    // This is necessary when we compute displacement field inplace as the scalar 
    // field becomes a vector field. It takes no time and won't hurt anyway :)
    // Note that we will temporarily store the density in 'potential' as the FFT will 
    // be computed in place ...
    potential.reinterpretNFields(); 
    potential.setName("density");    

    // cleanup previous data
    potential.erase();
    // Scatter projected density from AMR
    potential.addAmrGrid(&localAmrDensity);
    elapsed = scatterDensityTimer->stop();
    dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);
    
    if (symmetry != dice::RGSV::NONE)
      {
	scatterDensityTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Enforcing '%s' symmetry on density field ... ",symmetryStr.c_str());
	potential.getLocalGrid()->applySymmetry(symmetry);
	elapsed = scatterDensityTimer->stop();
	dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);
      }

    if (fileDumps.checkEvent(FileDumps::Density,true))
      {
	dumpTimer->start();
	potential.toVtk(fileDumps.getGlobalFName(FileDumps::Density).c_str(),
			fileDumps.getLocalFNameFormat(FileDumps::Density).c_str());
	double elapsed = dumpTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>
	  ("File was dumped in %lgs.\n",elapsed);
      }
    
    dumpGridDensityProfile();

    if (dice::glb::debug>2)
      {
	LocalGrid gatheredDensity;
	potential.gatherAll(gatheredDensity);
	gatheredDensity.
	  toVtk((fileDumps.getLocalFName(FileDumps::Potential)+
		 std::string("_gathered")).c_str());
      }    

    // save the projected density for computing potential energy later on (density is 
    // stored in 'potential' for now but will be transformed in place)
    typename dice::TimerPool::Timer timer;
    timer.start();
    dice::glb::console->printFlush<dice::LOG_INFO>("Cloning density ... ");
    clonedDensity.copyRawBuffer(potential);    
    dice::glb::console->printFlush<dice::LOG_INFO>("done in %lgs.\n",timer.stop());

    // Compute max density and setup density field if necessary
    if (units.useCosmo)
      {
	double volume=1.0L;
	for (int i=0;i<NDIM;++i)
	  volume *= mesh->getParams().delta[i];
	double avg = mass/volume;
	double fact = (3.0*a*units.cosmology.getParams().oM)/2.0;

	densityMax = (potential.getMax()/avg)*fact;	
	//FIXME: that's ugly, implement with functors !
	potential.subMultiply(avg,fact/avg);
      }
    else
      {
	densityMax = potential.getMax();
      }       

    // Finally, compute the potential or displacement field via FFT
    fftSolverTimer->start();

#ifdef D_FFT_POISSON_SOLVER_COMPUTE_DISPLACEMENT
    dice::glb::console->printFlush<dice::LOG_STD>("Computing displacement field ... ");
    fftSolver.execute();
    potential.setName("disp");
#else
    dice::glb::console->printFlush<dice::LOG_STD>("Computing potential ... ");
    fftSolver.execute();
    potential.setName("potential");    
#endif

    elapsed = fftSolverTimer->stop();
    dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);

    dice::glb::console->printFlush<dice::LOG_STD>("Computing grid statistics ... ");
    statisticsTimer->start();
          
    double oldPotentialEnergy=potentialEnergy;
    
#ifndef D_FFT_POISSON_SOLVER_COMPUTE_DISPLACEMENT
    //Compute potential energy
    PotentialEnergyGQT<RegularGrid,dice::gridKernel::InterpolateT<D,0>,0 >
      potentielEnergyGF(&clonedDensity,&potential);
    
    potentialEnergy = potential.quadrature(potentielEnergyGF)*unitsP;
#else
    dice::glb::console->printFlush<dice::LOG_STD>("(skipped Ep) ");
    potentialEnergy=0;    
#endif
    

    // Compute projected mass
    MassGQT< RegularGrid, dice::gridKernel::InterpolateT<D,0> >
      massGF(&clonedDensity);
    
    projectedMassErr = clonedDensity.quadrature(massGF);
    
    elapsed = statisticsTimer->stop();
    dice::glb::console->print<dice::LOG_STD>("done in %lgs.\n",elapsed);

    dice::glb::console->indent();	
    dice::glb::console->print<dice::LOG_INFO>
      ("-> Potential energy Ep = %e\n",potentialEnergy);
   
    dice::glb::console->print<dice::LOG_INFO>
      ("-> Projected mass pMass = %e (~%2.2e%% error)\n",
       projectedMassErr,100.0*std::abs(projectedMassErr-mass)/mass);
    projectedMassErr=std::abs(projectedMassErr-mass)/mass;	
    
    dice::glb::console->unIndent();

    if (fileDumps.checkEvent(FileDumps::Potential,true))
      {
	dumpTimer->start();
	potential.toVtk(fileDumps.getGlobalFName(FileDumps::Potential).c_str(),
			fileDumps.getLocalFNameFormat(FileDumps::Potential).c_str());

	double elapsed = dumpTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>
	  ("File was dumped in %lgs.\n",elapsed);
      }

    if (units.useCosmo)
      {
	double dt=0.5*(curSolverDeltaT+oldSolverDeltaT);
	double a0=units.cosmology.a_of_tau(t-dt,aStart);	
	double amid=units.cosmology.a_of_tau(t-0.5*curSolverDeltaT,aStart);

	expansionEnergy-=(a-a0)/amid *
	  (oldPotentialEnergy*oldSolverDeltaT+potentialEnergy*curSolverDeltaT)/
	  (oldSolverDeltaT+curSolverDeltaT);
      }
    else expansionEnergy=0;
  
    elapsed=solvePoissonTimer->stop();
    dice::glb::console->unIndent();
    dice::glb::console->print<dice::LOG_STD>("Done in %lgs.\n",elapsed);

    counter++;
  }
  
  template <class CT=IFloat>
  void updateProjectedDensityField()
  {
    typedef typename Vertex::ProjectedDensity::Type PDType;
    typedef typename Simplex::ProjectedDensityGradient::Type PDGType;
    int nThreads = dice::glb::num_omp_threads;

    // This is important if CT is doublebouble or quaddouble as FPU needs
    // to be set in the correct rounding mode for the results to be correct
    dice::hlp::FpuRoundingModeGuardT<CT> fpuGuard;
  
#ifdef D_TRADE_MEMORY_FOR_SPEED
    std::vector<double> vertexVolume(mesh->getNVertices()*2*nThreads,0);    
#else
    std::vector<double> vertexVolume(mesh->getNVertices(),0);
    
    // Reset vertices density
    FOREACH_THREAD_VERTEX(mesh,nThreads,th,it)
      {
	for (;it!=it_end;++it)
	  {
	    Vertex *v=(*it);
	    v->projectedDensity.setValue(0);
	  }
      }
    FOREACH_THREAD_END;
#endif
    
    // Compute the projected volume of each simplex
    FOREACH_THREAD_SIMPLEX_LG(mesh,nThreads,th,it)
      {
#ifdef D_TRADE_MEMORY_FOR_SPEED
	double *vVol = &vertexVolume[2*th];
	const long fac = 2*nThreads;
#endif
        for (;it!=it_end;++it)
	  {	    
	    Simplex *s=(*it);
	    //PDType volume = s->cache.d;
	    PDType mass = s->mass.getValue();

	    PDType volume = dice::hlp::numericStaticCast<PDType>
	      (mesh->template computeProjectedVolume<Simplex,CT>(s));
	    	    
	    /*
	    // This version interpolates the metric at vertices 
	    // KEEP COMMENTED FOR NOW
	    // This would only be usefull if we could project the tesselation of the
	    // quadratic reconstruction of the mesh instead of its raw linear version
	    
	    double volume[Simplex::NVERT];
	    dice::QuadraticSimplexT<NDIM,NDIM,double> 
	      qSimplex(s,s->segTracers.getConstPointer(),geometry);
	    	    
	    qSimplex.evalJacobianDetAtVertices(volume);
	    for (int i=0;i<Simplex::NVERT;++i) 
	      volume[i]=fabs(volume[i]);
	    */
	    
	    
	    for (int i=0;i<Simplex::NVERT;++i)
	      {	   
		Vertex *v=s->getVertex(i);
		if (v->isLocal())
		  {

#ifdef D_TRADE_MEMORY_FOR_SPEED

		    PDType *d=&vVol[v->getLocalIndex()*fac];
		    d[0]+=mass;
		    d[1]+=volume;

#else

		    PDType *tmp=v->projectedDensity.getPointer();
#pragma omp atomic update
		    (*tmp)+=mass;		    
#pragma omp atomic update
		    vertexVolume[v->getLocalIndex()]+=volume;

#endif
		  }
	      }
	    
	  }
      }
    FOREACH_THREAD_END;

#ifdef D_TRADE_MEMORY_FOR_SPEED

#pragma omp parallel num_threads(nThreads)
    {
      int th = omp_get_thread_num();
      const long fac = 2*nThreads;
      unsigned long start=fac*(mesh->getNVertices()/nThreads)*th;
      unsigned long stop =fac*(mesh->getNVertices()/nThreads)*(th+1);
      
      if (th==nThreads-1) stop=mesh->getNVertices()*fac;
      
      for (int i=start;i<stop;i+=fac)
	{
	  for (int j=2;j<fac;j+=2)
	    {
	      vertexVolume[i]+=vertexVolume[i+j];
	      vertexVolume[i+1]+=vertexVolume[i+j+1];
	    }
	}
    }

#endif

    // compute the density at each vertex
    FOREACH_THREAD_VERTEX(mesh,nThreads,th,it)
      {
#ifdef D_TRADE_MEMORY_FOR_SPEED
	const long fac = 2*nThreads;
#endif
	for(;it!=it_end;++it)
	  {
	    Vertex *v=(*it);
#ifdef D_TRADE_MEMORY_FOR_SPEED
	    PDType *d=&vertexVolume[v->getLocalIndex()*fac];
	    (*v->projectedDensity.getPointer()) = d[0]/d[1];
#else
	    PDType *d = v->projectedDensity.getPointer();
	    (*d)/=vertexVolume[v->getLocalIndex()];
#endif

	  }
      }
    FOREACH_THREAD_END;
  
    // and compute the gradient of the density for each simplex if necessary
    FOREACH_THREAD_SIMPLEX(mesh,nThreads,th,it)
      {
        for (;it!=it_end;++it)
	  {	    
            Simplex *s=(*it);
            PDType d[Simplex::NVERT];
            const Coord *vCoords[Simplex::NVERT];
	    for (int i=0;i<Simplex::NVERT;++i)
	      {
                Vertex *v=s->getVertex(i);
		d[i] = v->projectedDensity.getValue();
		vCoords[i]=v->getCoordsConstPtr();
	      }

	    PDGType *gradient=s->projectedDensityGradient.getPointer();

	    dice::BarycentricCoordinatesT<D,CT,GeometricProperties> 
	      bc(vCoords,*geometry);

	    bc.computeGradient(d,gradient);	
	  }
      }
    FOREACH_THREAD_END;
  }

  double getVelocityMax()
  {
    double vMax=0;
    std::vector<double> tmpMax(dice::glb::num_omp_threads,0);
    
    FOREACH_THREAD_VERTEX(mesh,dice::glb::num_omp_threads,th,it)
      {
	for (;it!=it_end;++it)
	  {
	    Vertex *v=(*it);
	    Coord *c=v->getCoordsPtr();
	    
	    for (int i=0;i<NDIM;++i) 
	      if (fabs(c[i+NDIM])>tmpMax[th]) tmpMax[th]=fabs(c[i+NDIM]);
	  }
      }
    FOREACH_THREAD_END;

    for (int j=0;j<dice::glb::num_omp_threads;j++) 
      if (tmpMax[j]>vMax) vMax = tmpMax[j];

    return vMax;
  }

  void updateInvariantThreshold(double t)
  {    
    double maxBoxLen=mesh->getParams().delta[0];
    for (int i=0;i<NDIM;++i)
      if (mesh->getParams().delta[i]>maxBoxLen)
	maxBoxLen=mesh->getParams().delta[i];

    double invariantFactor=units.velocity*units.length;
    if (units.useCosmo)
      {	
	double a=units.cosmology.a_of_tau(t,aStart);
	
	invariantFactor *= 
	  a/(units.H*maxBoxLen*maxBoxLen*units.length*units.length);	 
      }

    invariantThreshold_inv = invariantFactor/invariantThreshold;

    dice::glb::console->print<dice::LOG_STD>
      ("Updated invariant Threshold (th=%g): %g\n",
       invariantThreshold,1.0/invariantThreshold_inv);
  }   

  void updateTimeStep(double t, double &dt)
  {
    double upperLimit=dtMax;
    double oldDt=dt;
    double velocityMax = getVelocityMax();
    
    double cflFactor=units.length/units.velocity;
    const double pi=4.0*atan(1.0);
    double densityFactor=4.0*pi*units.G;
    double a=1;

    if (units.useCosmo)
      {
	a=units.cosmology.a_of_tau(t,aStart);
	double volume=1.0L;
	for (int i=0;i<NDIM;++i)
	  volume *= mesh->getParams().delta[i];
	double avgDensity = mass/volume;

	upperLimit=dtMax*a*units.cosmology.dtau_da(a); //dtMax = max(da/a) for cosmo
	cflFactor *= units.H/a;
	densityFactor=(3.0*a*units.cosmology.getParams().oM)/(2.0*avgDensity);
      }

    if (cflCondition == CflConditionTypeV::CFL_RHOMAX)
      {
	double dt1 = (cflRhoMax/sqrt(densityMax*densityFactor));
	double dt2 = (cflSizeMax/velocityMax)*cflFactor;
	
	dt = std::min(dt1,dt2);
	dt = std::min(dt,upperLimit);

	dt=mpiCom->min(dt);    
	dice::glb::console->print<dice::LOG_STD>
	  ("Updated time step (dtMax=%g): dt=%g->%g (rhoMax=%e,dt=%g)/(vMax=%e,dt=%g)\n",
	   upperLimit,oldDt,dt,densityMax,dt1,velocityMax,dt2);
      }
    else if (cflCondition == CflConditionTypeV::RHOMAX)
      {
	dt = (cflRhoMax / sqrt(densityMax*densityFactor));
	dt = std::min(dt,upperLimit);

	dt=mpiCom->min(dt);    
	dice::glb::console->print<dice::LOG_STD>
	  ("Updated time step (dtMax=%g): dt=%g->%g (rhoMax=%e)\n",
	   upperLimit,oldDt,dt,densityMax);
      }
    else if (cflCondition == CflConditionTypeV::CFL)
      {
	if (velocityMax==0) 
	  dt = upperLimit;
	else
	  dt = (cflSizeMax/velocityMax)*cflFactor;

	dt = std::min(dt,upperLimit);

	dt=mpiCom->min(dt);    
	dice::glb::console->print<dice::LOG_STD>
	  ("Updated time step (dtMax=%g): dt=%g->%g (vMax=%e)\n",
	   upperLimit,oldDt,dt,velocityMax);
      }   

    if (units.useCosmo)
      {
	dice::glb::console->print<dice::LOG_STD>
	  ("New expansion factor is a=%g / z=%g (da=%g)\n",
	   a,1.0/a-1.0,units.cosmology.a_of_tau(t+dt,aStart)-a);
      }
  }

  void drift_slow(double dt)
  {
    // Multithreading is limited by memory bandwidth here and having too many thread
    // may result in a significant slowdown so we need to put a limit
    int nThreads=std::min(dice::glb::num_omp_threads,D_N_THREADS_MAX_DRIFT);
    double dth=dt/2;

    //MeshAndTracersCoords coordContainer(mesh);
    //int nThreads = dice::glb::num_omp_threads;
    MeshAndTracersCoords coordContainer(mesh);
    //int nBatches = nThreads;
#pragma omp parallel num_threads(nThreads) 
      {	
	int th = omp_get_thread_num();
	const MeshAndTracersCoords_iterator it_end = 
	  coordContainer.end(th,nThreads);

	for (MeshAndTracersCoords_iterator it = 
	       coordContainer.begin(th,nThreads);it!=it_end;++it)
	  {	   
	    Coord *c=*it;
	    
	    for (int i=0;i<NDIM;++i)
	      c[i]+=c[i+NDIM]*dth;
	    
	    geometry->sanitizeBoundary(c);
	  }
      }    
  }

  void drift(double dt)
  {
    // Multithreading is limited by memory bandwidth here and having too many thread
    // may result in a significant slowdown so we need to put a limit
    int nThreads=std::min(dice::glb::num_omp_threads,D_N_THREADS_MAX_DRIFT);    

#pragma omp parallel num_threads(nThreads)
    {LIKWID_MARKER_START("Drift");}
        
    double dth=dt/2;
    typename dice::TimerPool::Timer timer;
    
    double driftFactor=1.0/(units.length/units.velocity);
    if (units.useCosmo)
      {
	double a=units.cosmology.a_of_tau(curSolverTime,aStart);
	// h0 -> H0 in s-1	
	driftFactor*=a/units.H;
      }
    else
      {
	double T=1.0;
	driftFactor*=T;	
      }
    
    timer.start();
#pragma omp parallel num_threads(nThreads) 
    {
      int th = omp_get_thread_num();
      const vertexPtr_LGS_iterator itv_end=mesh->vertexLGSEnd();
      for (vertexPtr_LGS_iterator it=mesh->vertexLGSBegin(th,nThreads);
	   it!=itv_end;++it)
	{
	  Coord *c=it->getCoordsPtr();
	  for (int i=0;i<NDIM;++i)
	    {
	      c[i]+=c[i+NDIM]*driftFactor*dth;	    
	    }
	  geometry->sanitizeBoundary(c);
	}
    }
    
#pragma omp parallel num_threads(nThreads) 
    {
      int th = omp_get_thread_num();
      const simplexPtr_LGS_iterator its_end=mesh->simplexLGSEnd();
      for (simplexPtr_LGS_iterator it=mesh->simplexLGSBegin(th,nThreads);
	   it!=its_end;++it)
	{
#ifndef NO_SIMPLEX_TRACERS	  	  
	  static const int nTracers = (Mesh::Simplex::NSEG+1);
#else
	  static const int nTracers = (Mesh::Simplex::NSEG);
#endif

	  DICE_ALIGNAS(32) Coord temp1[nTracers*NDIM];
	  DICE_ALIGNAS(32) Coord temp2[nTracers*NDIM];
	  int k=0;

#ifndef NO_SIMPLEX_TRACERS
	  Coord *c=it->tracer.getPointer();
	  for (int j=0;j<NDIM;++j,++k)
	    {
	      temp1[k]=c[j];
	      temp2[k]=c[j+NDIM];
	    }
#endif	  

	  Coord *c2=it->segTracers.getPointer();	  	  
	  for (int i=0;i<Mesh::Simplex::NSEG*NDIM_W;i+=NDIM_W)
	    for (int j=0;j<NDIM;++j,++k)
	      {
		temp1[k]=c2[i+j];
		temp2[k]=c2[i+j+NDIM];
	      }	    

	  for (int i=0;i<nTracers*NDIM;++i)
	    temp1[i]+=temp2[i]*driftFactor*dth;
	  for (int i=0;i<nTracers*NDIM;i+=NDIM)
	    {	      
	      geometry->sanitizeBoundary(temp1+i);
	    }
	    
	  k=0;
#ifndef NO_SIMPLEX_TRACERS
	  for (int j=0;j<NDIM;++j,++k)
	    c[j]=temp1[k];
#endif	    
	  for (int i=0;i<Mesh::Simplex::NSEG*NDIM_W;i+=NDIM_W)
	    for (int j=0;j<NDIM;++j,++k)
	      c2[i+j]=temp1[k];
	    
	}
    }
    

#pragma omp parallel num_threads(nThreads)
    {LIKWID_MARKER_STOP("Drift");}
  }

  void kick_drift(double dt)
  {
    MeshAndTracersCoords cc(mesh);
    kick_drift(dt,cc);
  }

  template <class CC>
  void kick_drift(double dt, CC &coordContainer)
  {
#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_START("KickDrift");}
   
    double dth=dt/2;
    //MeshAndTracersCoords coordContainer(mesh);
    
    double dispFactor=(units.length/units.velocity);
    double velFactor=1.0;
    double driftFactor=1.0;

    if (units.useCosmo)
      {
	double a0=units.cosmology.a_of_tau(curSolverTime,aStart);
	double a1=units.cosmology.a_of_tau(curSolverTime + curSolverDeltaT,aStart);
	 
	velFactor = a0/a1;
	// h0 -> H0 in s-1	
	dispFactor *= units.H/a1;
	driftFactor = 1.0/dispFactor;
      }
    else
      {
	double T=1.0;
	dispFactor/=T;
	driftFactor = 1.0/dispFactor;
      }
    
    
#pragma omp parallel for num_threads(dice::glb::num_omp_threads)
    for (int j=0;j<dice::glb::num_omp_threads;j++)
      {
     	double disp[NDIM];
	InterpolationKernel interpolator;
	
	gatheredPotential.initializeInterpolationKernel(interpolator);        

	const auto it_end = coordContainer.end(j,dice::glb::num_omp_threads);
	for (auto it = coordContainer.begin(j,dice::glb::num_omp_threads);it!=it_end;++it)
	  {	    
	    Coord *c=*it;
	    gatheredPotential.applyKernel(interpolator,c,disp,-1);

	    for (int i=0;i<NDIM;++i)
	      c[i+NDIM]=c[i+NDIM]*velFactor + disp[i]*dispFactor*dt;
	    for (int i=0;i<NDIM;++i)
	      c[i]+=c[i+NDIM]*driftFactor*dth;

	    geometry->sanitizeBoundary(c);
	  }
      } 

#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_STOP("KickDrift");}
  }

  void advance_impl()
  { 
    double elapsed;

    if (dice::glb::debug>2)
      mesh->dumpToNDnetwork((fileDumps.getLocalFName(FileDumps::Mesh)+
			     std::string("_beforeAdvance")).c_str(),
			    dice::IO::NDNET_WithShadows|
			    dice::IO::NDNET_WithGhosts|
			    dice::IO::NDNET_WithNeighbors,
			    true);
    
    dice::glb::console->printFlush<dice::LOG_STD>("Advancing (drift) ... ");
    driftTimer->start();
    drift(curSolverDeltaT);
    elapsed=driftTimer->stop();
    dice::glb::console->printFlush<dice::LOG_STD>("done in %.3gs.\n",elapsed);

    if (dice::glb::debug>2)
      mesh->dumpToNDnetwork((fileDumps.getLocalFName(FileDumps::Mesh)+
			     std::string("_afterDrift")).c_str());
        
    solvePoisson(curSolverTime + curSolverDeltaT/2);

#ifdef D_FFT_POISSON_SOLVER_COMPUTE_DISPLACEMENT   
    dice::glb::console->printFlush<dice::LOG_STD>("Gathering displacement field ... ");
#else
    dice::glb::console->printFlush<dice::LOG_STD>("Gathering potential ... ");
#endif
    
    
    double maxAlloc=gatheredPotentialAllocLimit*(1L<<30);
    if (MeshAndTracersSliceCoordsT<Mesh>::getSlicesCount(potential,maxAlloc)>1)
      {	
	// We cannot allocate a full grid locally, so we need several passes
	MeshAndTracersSliceCoordsT<Mesh> coordContainer(mesh,potential,maxAlloc);
	dice::glb::console->print<dice::LOG_STD>("\n");
	dice::glb::console->indent();
	int pass=0;
	int nRequired=0;
	do {
	  coordContainer.setGroup(pass);
	  dice::glb::console->printFlush<dice::LOG_STD>("Pass %d: gathering ...",pass++);
	  gatherPotentialTimer->start();
	  potential.template gatherSubsetAtCoords<InterpolationKernel>
	    (coordContainer,gatheredPotential);
	  double elapsedG = gatherPotentialTimer->stop();
	  dice::glb::console->printFlush<dice::LOG_STD>(" advancing (K+D) ... ");
	  kickAndDriftTimer->start();
	  kick_drift(curSolverDeltaT);
	  double elapsedK=kickAndDriftTimer->stop();

	  nRequired = 
	    mpiCom->sum(std::distance(coordContainer.begin(),coordContainer.end())!=0);

	  dice::glb::console->print<dice::LOG_STD>
	    ("done in %.3gs / %.3gs. (%d/%d required)\n",
	     elapsedG,elapsedK,nRequired,mpiCom->size());
	} while ((nRequired!=0)&&(pass<coordContainer.nGroups()));
	dice::glb::console->unIndent();
      }
    else
      {
	// We can afford to allocate the full grid locally
	gatherPotentialTimer->start();
    	potential.template gatherSubsetAtCoords<InterpolationKernel>
	  (MeshAndTracersCoords(mesh),gatheredPotential);    
 
	elapsed = gatherPotentialTimer->stop();
	dice::glb::console->printFlush<dice::LOG_STD>("done in %lgs.\n",elapsed);  

	if ((dice::glb::debug)&&(fileDumps.recheckEvent(FileDumps::Potential)))
	  {
	    gatheredPotential.
	      toVtk((fileDumps.getLocalFName(FileDumps::Potential)+
		     std::string("_atCoords")).c_str());
	  }    

	dice::glb::console->printFlush<dice::LOG_STD>("Advancing (kick+drift) ... ");
	kickAndDriftTimer->start();
	kick_drift(curSolverDeltaT);
	elapsed=kickAndDriftTimer->stop();
	dice::glb::console->printFlush<dice::LOG_STD>("done in %.3gs.\n",elapsed);
      }
    
    dumpNetwork();   
  }  

protected:
  
  void dumpGridDensityProfile(bool force=false)
  {
    if ((force)||(fileDumps.checkEvent(FileDumps::RadialGridDensity,true)))
      {
	statisticsTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Computing grid density profile ... ");

	RadialProfileGT<RegularGrid> profileG(&potential);
	if (specifyProfileCenter) profileG.setCenterCoordinates(densityProfileCenter);
	potential.visit(profileG);
	double elapsed=statisticsTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>("done in %lgs.\n",elapsed);

	dumpTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Dumping grid density profile to file '%s' ... ",
	   fileDumps.getGlobalFName(FileDumps::RadialGridDensity).c_str());

	char comment[255];
	sprintf(comment,"t=%lg step=%ld",curSolverTime,curStepIndex);
	profileG.
	  toFile(fileDumps.getGlobalFName(FileDumps::RadialGridDensity),mpiCom,comment);
	elapsed = dumpTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>("done.\n");
      }
  }

  void dumpMeshDensityProfile(bool force=false, int nSamplesPerOvelap=10)
  {
    if ((force)||(fileDumps.checkEvent(FileDumps::RadialMeshDensity,true)))
      {
	statisticsTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Computing mesh density profile ... ");

	// Compute r0 as 1/5th of a pixel extent
	double x0[NDIM];
	double delta[NDIM];
	mesh->getBoundingBox(x0,delta,false);
	double r0=delta[0]/(1<<fftGridLevel);
	for (int i=1;i<NDIM;++i)
	  {
	    double tmp=delta[i]/(1<<fftGridLevel);
	    if (r0<tmp) r0=tmp;
	  }
	r0*=0.2;       	
	
	RadialProfileMT<Mesh> profileM(mesh,nSamplesPerOvelap,r0,true);
	if (specifyProfileCenter) profileM.setCenterCoordinates(densityProfileCenter);
	mesh->visitSimplices(profileM,true,false,false);
	double elapsed=statisticsTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>("done in %lgs.\n",elapsed);

	dumpTimer->start();
	dice::glb::console->printFlush<dice::LOG_STD>
	  ("Dumping mesh density profile to file '%s' ... ",
	   fileDumps.getGlobalFName(FileDumps::RadialMeshDensity).c_str());

	char comment[255];
	sprintf(comment,"t=%lg step=%ld",curSolverTime,curStepIndex);
	profileM.
	  toFile(fileDumps.getGlobalFName(FileDumps::RadialMeshDensity),mpiCom,comment);
	elapsed = dumpTimer->stop();
	dice::glb::console->printFlush<dice::LOG_INFO>("done.\n");
      }
  }

  void dumpNetwork(bool forceFull=false,bool forceSubNetwork=false)
  {
    dumpTimer->start();
    int dumped=false;
	
    if ((forceSubNetwork)||(fileDumps.checkEvent(FileDumps::Caustics,true)))
      {dumped=true;dumpCaustics();}

    if ((forceSubNetwork)||(fileDumps.checkEvent(FileDumps::Subsets,true)))
      {dumped=true;dumpSubsets();}

    if ((forceSubNetwork)||(fileDumps.checkEvent(FileDumps::Lines,true)))
      {dumped=true;dumpLagrangianLines();}

    if ((forceFull)||(fileDumps.checkEvent(FileDumps::Mesh,true)))
      {dumped=true;mesh->dumpToNDnetwork(fileDumps.getLocalFName(FileDumps::Mesh).c_str());}

    bool dumpedFatMesh=false;
    if ((forceFull)||(fileDumps.checkEvent(FileDumps::FatMesh,true)))
      {
	dumped=true;
	auto *fptr=mesh->getSimplexFunctorPtr("segTracers");
	long flags=0;

	if (fptr!=NULL)
	  {
	    // Set flags to enable dumping segment tracers
	    flags=fptr->getFlags();
	    fptr->setFlags(dice::cellDataFunctors::F_NO_FLAG);
	  }
	
	mesh->dumpToNDnetwork(fileDumps.getLocalFName(FileDumps::Mesh).c_str(),
			      dice::IO::NDNET_WithShadows|
			      dice::IO::NDNET_WithGhosts|
			      dice::IO::NDNET_WithNeighbors);
	
	if (fptr!=NULL) fptr->setFlags(flags);
	
	dumpedFatMesh=true;
      }
    
    if ((forceFull)||(fileDumps.checkEvent(FileDumps::Mesh,true)))
      {
	if (!dumpedFatMesh)
	  {
	    dumped=true;
	    mesh->dumpToNDnetwork(fileDumps.getLocalFName(FileDumps::Mesh).c_str());
	  }
      }
          
    double elapsed = dumpTimer->stop();

    if (dumped) dice::glb::console->printFlush<dice::LOG_INFO>
		  ("Files were dumped in %lgs.\n",elapsed);
  }
   
  void dumpCaustics()
  {
    auto causticsFilter = [&](const FacetHandle &f)->FacetHandle 
      {
	if ((mesh->computeProjectedVolume(f->getSimplex(),true)*
	     mesh->computeProjectedVolume(f->getOppositeSimplex(),true))<0)
	  return f;
	else return FacetHandle();
      };

    mesh->template dumpToFilteredNDNetwork<decltype(causticsFilter),FacetHandle>
      (fileDumps.getLocalFName(FileDumps::Caustics).c_str(),
       causticsFilter,dice::IO::NDNET_NoSimplices|dice::IO::NDNET_TrimExtraVertices);
  }

  void dumpSubsets()
  {
    const char *patterns[]={"*","_*_","__*__","___*___","____*____"};
    const char *p[NDIM];

    if (NDIM==2)
      for (int i=0;i<NDIM;++i)
	p[i]=patterns[2*i];
    else
      for (int i=0;i<NDIM;++i)
	p[i]=patterns[NDIM-1];
    
    NDnetFilter_subsetsT<Mesh> subsetsFilter(mesh,initialMeshResolution,p);
   
    mesh->template dumpToFilteredNDNetwork<decltype(subsetsFilter),Simplex>
      (fileDumps.getLocalFName(FileDumps::Subsets).c_str(),
       subsetsFilter,dice::IO::NDNET_TrimExtraVertices);
  }

  void dumpLagrangianLines()
  {     
    typedef NDnetFilter_lagrangianLinesT<Mesh> LinesFilter;
    typedef typename LinesFilter::Handle LinesFilterHandle;

    int delta[]={2*NDIM-1,2*NDIM-1,2*NDIM-1};
    LinesFilter linesFilter(mesh,initialMeshResolution,delta);

    mesh->template dumpToFilteredNDNetwork<LinesFilter,LinesFilterHandle>
      (fileDumps.getLocalFName(FileDumps::Lines).c_str(),
       linesFilter,dice::IO::NDNET_NoSimplices|dice::IO::NDNET_TrimExtraVertices);   
  }

  // All the variables we use are down there !
protected:
  IC *ic;
  Mesh *mesh;  
  GeometricProperties *geometry;
  LocalAmrGrid localAmrDensity;  
  RegularGrid clonedDensity;
  RegularGrid potential;
  LocalGrid gatheredPotential;
  FFTSolver fftSolver;  

  dice::MpiCommunication *mpiCom;

  typename dice::TimerPool::Timer *amrBuildTimer;
  typename dice::TimerPool::Timer *projectTimer;   
  typename dice::TimerPool::Timer *projectBarrierTimer;  
  typename dice::TimerPool::Timer *kickAndDriftTimer;  
  typename dice::TimerPool::Timer *driftTimer;  
  typename dice::TimerPool::Timer *fftSolverTimer;
  typename dice::TimerPool::Timer *updateDensityTimer;
  typename dice::TimerPool::Timer *solvePoissonTimer;  
  typename dice::TimerPool::Timer *scatterDensityTimer;
  typename dice::TimerPool::Timer *gatherPotentialTimer;
  typename dice::TimerPool::Timer *sortMeshTimer;
  typename dice::TimerPool::Timer *statisticsTimer;
  typename dice::TimerPool::Timer *dumpTimer;

  long curStepIndex;
  double curSolverTime;
  double curSolverDeltaT;
  double oldSolverDeltaT;
  
  double aStart;
  double aEnd;
  double dtMax;

  double mass;  
  Units units;
  
  float serializedVersion;

  CflConditionType cflCondition;  
  double cflGrid;
  double cflRhoMax;
  double cflSizeMax;

  dice::RegularGridSymmetryE symmetry;
  std::string symmetryStr;

  int initialMeshResolution[NDIM];

  double refineThreshold;
  double coarsenHysteresis;
  double volumeThreshold;
  double invariantThreshold;

  double refineThreshold2;
  double coarsenThreshold2; 

  double volumeThreshold_inv;
  double invariantThreshold_inv;
  double maximumSegmentLength2_inv;  

  double phSortThreshold;

  int checkProjectedDensity;
  double accuracyLevel;

  double gatheredPotentialAllocLimit;
  int fftGridLevel;
  int fftWisdom;
  std::string importFftWisdom;
  std::string exportFftWisdom;
  
  int squeeze;
  int skipInitialPoisson;
  int dumpInitialMesh;
  int noRepartWeight;
  int projectionOrder;

  int splitLongestEdge;
  int maxAmrLevel;
  int fastAmrBuild;

  int pForceHighResolution;
  int pResetHRMode;
  //std::string pSwitchHRModeFileName;
  double pEnableHRModeThreshold;
  double pAnisotropyThreshold;
  double pWidthThreshold;

  int pUseVerticesThreshold;
  int maxSimplexLevel;
  
  int specifyProfileCenter;
  double densityProfileCenter[NDIM];
 
  long newMeshSimplicesCount; 
  
  double densityMax;

  double projectedMassErr;
  double totalVolume;
  double totalVolume_Q;
  double kineticEnergy;
  double kineticEnergy_Q;
  double potentialEnergy; 

  double expansionEnergy;
 
  bool repartStatus;
  int rebuildAmrEvery;

  FileDumps fileDumps;
};

#endif
