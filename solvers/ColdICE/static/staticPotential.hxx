#ifndef __COLDICE_STATIC_POTENTIAL_HXX__
#define __COLDICE_STATIC_POTENTIAL_HXX__

#include <dice/cosmo/staticPotentialSolutions.hxx>

// This is just a convient hack to compute statistics from a snapshot. 

template <class VlasovPoissonSolver>
class StaticPotential : public VlasovPoissonSolver
{
public:
  typedef VlasovPoissonSolver Base;
  typedef StaticPotential<VlasovPoissonSolver> MyType;

  typedef dice::cosmo::StaticPotentialSolutionsT<Base::NDIM> StaticPotentialSolutions;
  typedef typename StaticPotentialSolutions::PotentialType  PotentialType;
  typedef typename StaticPotentialSolutions::PotentialTypeV PotentialTypeV;
  typedef typename StaticPotentialSolutions::PotentialTypeSelect PotentialTypeSelect;

  typedef typename Base::Mesh Mesh;
  typedef typename Base::MeshParams MeshParams;

  typedef typename Base::SimplexFunctor SimplexFunctor;
  typedef typename Base::VertexFunctor  VertexFunctor;

  typedef typename Base::Coord Coord;
  typedef typename Base::MeshAndTracersCoords MeshAndTracersCoords;
  typedef typename Base::MeshAndTracersCoords_iterator MeshAndTracersCoords_iterator;
 
  static const int NDIM          = Base::NDIM;
  static const int NDIM_W        = Base::NDIM_W;
  static const int BOUNDARY_TYPE = Base::BOUNDARY_TYPE;

  static std::string parserCategory() {return "staticSolver";}
  // static std::string parserCategory() {return "post";}
  // static std::string classHeader() {return "vlasov_poisson_post";}
  // static float classVersion() {return 0.10;}
  // static float compatibleSinceClassVersion() {return 0.10;}

  template <class SP, class R, class PM>
  StaticPotential(MeshParams &meshParams,
		  SP &solverInterfaceParams,
		  R *reader,
		  PM& paramsManager,
		  float serializedVersion,
		  dice::MpiCommunication *mpiCom_):
    Base(meshParams,solverInterfaceParams,reader,
	 paramsManager,serializedVersion,mpiCom_)
  {   
    //typename PM::Parser *parser = paramsManager.getParser();

    potentialType = PotentialTypeV::PLUMMER;
    std::string potentialTypeStr = paramsManager.
      get("potential",parserCategory()
	  ,PotentialTypeSelect().getString(potentialType,true),reader,
	  PM::PARSER_FIRST,
	  PotentialTypeSelect().getAllString("Static potential type (%s)"));
    potentialType = PotentialTypeSelect().getVal(potentialTypeStr,true);

    gravitationalMass = 1.0;
    gravitationalMass = paramsManager.
      get("gravitationalMass",parserCategory(),gravitationalMass,reader,
	  PM::FILE_FIRST,
	  "Total mass generating the static potential");   
    /*
    radius = 0.2;
    radius = paramsManager.
      get("radius",parserCategory(),radius,reader,
	  PM::FILE_FIRST,
	  "Effective radius of the static potential");
    */

    analyticSolution.selectPotential(potentialType,reader,paramsManager,
				     gravitationalMass,Base::units.G);
  }
  /*
  template <class SP>
  void onInitialize(Mesh *m, const SP& solverParams, double t, bool restart)
  { 
    init_impl(m,solverParams,t,restart);
    Base::init_impl(m,solverParams,t,restart);     
  }
  */
  void onAdvance()
  {       
    advance_impl();
  }
 
protected:
  /*
  template <class SP>
  void init_impl(Mesh *m, const SP& solverParams, double t, bool restart)
  {  
    analyticSolution.setPotential(potentialType,gravitationalMass,Base::units.G);
  }  
  */

  void drift(double dt)
  {
#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_START("Drift");}

    double dth=dt/2; 
    
    MeshAndTracersCoords coordContainer(Base::mesh);
   
#pragma omp parallel for num_threads(dice::glb::num_omp_threads)
    for (int j=0;j<dice::glb::num_omp_threads;j++)
      {
   	const MeshAndTracersCoords_iterator it_end = 
	  coordContainer.end(j,dice::glb::num_omp_threads);

	for (MeshAndTracersCoords_iterator it = 
	       coordContainer.begin(j,dice::glb::num_omp_threads);it!=it_end;++it)
	  {	    
	    Coord *c=*it;
	    
	    for (int i=0;i<NDIM;++i)
	      c[i]+=c[i+NDIM]*dth;

	    Base::geometry->checkBoundary(c);
	  }
      }	
    
#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_STOP("Drift");}
  }

  void kick_drift(double dt)
  {
#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_START("KickDrift");}

    double dth=dt/2; 
    
    MeshAndTracersCoords coordContainer(Base::mesh);
    // Drift + Kick
#pragma omp parallel for num_threads(dice::glb::num_omp_threads)
    for (int j=0;j<dice::glb::num_omp_threads;j++)
      {
   	const MeshAndTracersCoords_iterator it_end = 
	  coordContainer.end(j,dice::glb::num_omp_threads);

	for (MeshAndTracersCoords_iterator it = 
	       coordContainer.begin(j,dice::glb::num_omp_threads);it!=it_end;++it)
	  {	    
	    Coord *c=*it;
	    double r=0;

	    for (int i=0;i<NDIM;++i) 
	      r+=c[i]*c[i];
	    r=sqrt(r);

	    //Drift
	    for (int i=0;i<NDIM;++i)
	      c[i+NDIM]  += analyticSolution.acc_over_r(i,r,c)*c[i]*dt; 

	    // Kick
	    for (int i=0;i<NDIM;++i)
	      c[i]+=c[i+NDIM]*dth;

	    Base::geometry->checkBoundary(c);
	  }
      } 

#pragma omp parallel num_threads(dice::glb::num_omp_threads)
    {LIKWID_MARKER_STOP("KickDrift");}
  }

  void advance_impl()
  {   
    bool computePoisson = 
      Base::fileDumps.checkEvent(FileDumps::Amr)||
      Base::fileDumps.checkEvent(FileDumps::Density)||
      Base::fileDumps.checkEvent(FileDumps::Potential);

    // Drift
    dice::glb::console->printFlush<dice::LOG_STD>("Advancing (D) ... ");
    Base::driftTimer->start();
    drift(Base::curSolverDeltaT);
    double elapsed=Base::driftTimer->stop();
    dice::glb::console->printFlush<dice::LOG_STD>("done in %.3gs.\n",elapsed); 

    // Solve poisson if we need the projected density or AMR grid ...
    if (computePoisson) Base::solvePoisson(Base::curSolverTime + Base::curSolverDeltaT/2);

    if ((dice::glb::debug)&&
	(Base::fileDumps.recheckEvent(FileDumps::Potential)))
      {	
	Base::potential.template gatherSubsetAtCoords<typename Base::InterpolationKernel>
	  (MeshAndTracersCoords(Base::mesh),Base::gatheredPotential);    
	Base::gatheredPotential.
	  toVtk((Base::fileDumps.getLocalFName(FileDumps::Potential)+
		 std::string("_atCoords")).c_str());
      } 

    if (Base::fileDumps.checkEvent(FileDumps::Mesh)||
	Base::fileDumps.checkEvent(FileDumps::Caustics)||
	Base::fileDumps.checkEvent(FileDumps::Lines)||
	Base::fileDumps.checkEvent(FileDumps::Subsets))
      {
	Base::updateProjectedDensityField();	
	Base::dumpNetwork();	
      }
    
    // Kick + Drift
    dice::glb::console->printFlush<dice::LOG_STD>("Advancing (KD) ... ");
    Base::kickAndDriftTimer->start();
    kick_drift(Base::curSolverDeltaT);
    elapsed=Base::kickAndDriftTimer->stop();
    dice::glb::console->printFlush<dice::LOG_STD>("done in %.3gs.\n",elapsed);    
       
  }
  
private:
  StaticPotentialSolutions analyticSolution;

  PotentialType potentialType;
  double gravitationalMass;
  //double radius;
};

#endif
