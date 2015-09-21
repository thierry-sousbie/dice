#ifndef __SOLVER_INTERFACE_HXX__
#define __SOLVER_INTERFACE_HXX__

#include "../dice_globals.hxx"

// This is the interface to the solver implementation 
#include <stdio.h>
#include <math.h>
#include <string>
#include <limits> 

#include "./solverInterfaceParams.hxx"

#include "../tools/helpers/helpers.hxx"
#include "../tools/IO/myIO.hxx"
#include "../tools/IO/IOHelpers.hxx"
#include "../tools/IO/paramsManager.hxx"
#include "../tools/wrappers/likwidWrapper.hxx"
#include "../tools/wrappers/pProfWrapper.hxx"

#include "../internal/namespace.header"

namespace slv {

template <class SI>
class SolverInterfaceT
{
private:

  enum {signal_STOP=(1<<0), 
	signal_DUMP_RESTART=(1<<1)};
 
public:

  typedef SolverInterfaceT<SI> MyType;

  typedef SI SolverImplementation;   
  typedef typename SI::Mesh Mesh;  

  typedef typename SolverImplementation::RefinementStatus RefinementStatus;
  typedef typename SolverImplementation::CoarseningStatus CoarseningStatus;

  typedef SolverInterfaceParamsT<RefinementStatus,CoarseningStatus> Params;

  typedef typename Mesh::Params MeshParams;
  typedef typename Mesh::Simplex Simplex;

  typedef myIO::BinaryReaderT<> BinaryReader;
  typedef ParamsManagerT<ParamsParser> ParamsManager;//,Console,LOG_INFO> ParamsManager;
 
  static std::string parserCategory() {return "solver";}
  static std::string classHeader() {return "generic_solver_base";}
  static float classVersion() {return 0.11;}
  static float compatibleSinceClassVersion() {return 0.11;}

  SolverInterfaceT(MpiCommunication *com=glb::mpiComWorld):
    manager(glb::pParser)//, glb::console),
    
  {    
    reader = BinaryReader::nullReader();
    mpiCom = com;
    glb::console->print<LOG_STD>("Initializing solver:\n");
    glb::console->indent();
    fromRestartFile=false;
    initRestart();
    construct();   
  }

  SolverInterfaceT(const std::string &restartFileName, 
		   MpiCommunication *com=glb::mpiComWorld):
    manager(glb::pParser),
    mpiCom(com)
  {    
    if (restartFileName.length() != 0)//std::string(""))
      {
	reader = new BinaryReader(formatRestartFileName(restartFileName));
	glb::console->print<LOG_STD>("Restarting solver from file '%s':\n",
				     formatRestartFileName(restartFileName).c_str());
	fromRestartFile=true;
      }
    else
      {
	reader = BinaryReader::nullReader(); 
	glb::console->print<LOG_STD>("Initializing solver:\n");
	fromRestartFile=false;
      }
    glb::console->print<LOG_STD>("\n");
    glb::console->indent();
    initRestart();
    construct();
  }

  ~SolverInterfaceT()
  {
    delete implementation;
    delete mesh;
    if (reader != NULL)
      delete reader;
  }

  void construct()
  {                   
    float implVersion;
    float paramsVersion;
    float meshParamsVersion;

    BinaryReader::checkHeaderAndReport<LOG_ERROR,LOG_WARNING,Params>
      (glb::console,reader,paramsVersion,true);
    params.setSerializedVersion(paramsVersion);

    BinaryReader::checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MeshParams>
      (glb::console,reader,meshParamsVersion,true);
    meshParams.setSerializedVersion(meshParamsVersion);

    BinaryReader::checkHeaderAndReport<LOG_ERROR,LOG_WARNING,SolverImplementation>
      (glb::console,reader,implVersion,true);

    implementation = new SolverImplementation
      (meshParams,params,reader,manager,implVersion,mpiCom);

    //setupParameters();
    setupNonParameters();

    //mesh = new Mesh(meshParams,reader,mpiCom,glb::pParser);
    mesh = new Mesh(mpiCom);

    globalTimer  = glb::timerPool->pop("global");
    buildTimer   = glb::timerPool->pop("build");
    stepTimer    = glb::timerPool->pop("step");
    advanceTimer = glb::timerPool->pop("advance");
    coarsenTimer = glb::timerPool->pop("coarsen");
    refineTimer  = glb::timerPool->pop("refine");
    repartTimer  = glb::timerPool->pop("repart");
    updateTimer  = glb::timerPool->pop("update");
    writeRestartTimer  = glb::timerPool->pop("writeRestart");

    //const MeshParams& p = mesh->getParams();

    glb::console->print<LOG_INFO>("\n");
    manager.reportParsed<LOG_INFO>();
    glb::console->print<LOG_INFO>("\n");

    glb::console->print<LOG_INFO>("The following managed parameters will be used :\n");
    glb::console->indent();
    manager.report<LOG_INFO>();
    glb::console->print<LOG_INFO>("\n");   
    glb::console->unIndent();   

    glb::console->print<LOG_INFO>("\n");
    
    if (curTime>=params.tEnd)
      {
	glb::console->print<LOG_WARNING>(" tEnd (=%lg) < curTime (=%lg)\n",params.tEnd, curTime);
	glb::console->print<LOG_WARNING>(" solver will return after making a single step !\n");
	params.tEnd = nextafter(curTime,std::numeric_limits<double>::max());
	//glb::console->print<LOG_WARNING>(" Are you sure that's what you want ??? You should probably change solver.tEnd ;)\n");
      }

    // This is for when we want to use likwid perfmon for profiling ...
    LIKWID_MARKER_INIT;
#pragma omp parallel num_threads(glb::num_omp_threads)
    {LIKWID_MARKER_THREADINIT;}  
  }

  std::string formatRestartFileName(const std::string &fname, bool write=false)
  {
    char buf[1024];
    if (write)
      {
	if (mpiCom->size()>1) 
	  sprintf(buf,"%s_%6.6d.rst",fname.c_str(),mpiCom->rank());
	else
	  sprintf(buf,"%s.rst",fname.c_str());
      }
    else
      {
	if (mpiCom->size()>1) 
	  sprintf(buf,"%s_%6.6d.rst",fname.c_str(),mpiCom->rank());
	else
	  sprintf(buf,"%s.rst",fname.c_str());
      }
    return std::string(buf);
  }

  void initRestart()
  {
    if (reader == NULL) 
      {
	version = classVersion();
	return;
      }
    
    BinaryReader::checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);

    if (reader->getNeedSwap())
      {
	glb::console->print<LOG_INFO>("Endianness difference detected, will swap binary data (slower).\n");
      }

    reader->checkTypeSize<LOG_WARNING_SINGLE>(glb::console);
    
  }

  void writeRestartFile(const std::string &fname, bool force=false)
  {    
    std::string fullFName(formatRestartFileName(params.outputDir+fname,true));    
    glb::console->printFlush<LOG_STD>("Writing restart file '%s' ... ",fullFName.c_str());
    
    if ((params.noRestart)&&(!force))
      {
	glb::console->printFlush<LOG_STD>("SKIPPED.\n");
	return;
      }

    writeRestartTimer->start();
    myIO::BinaryWriterT<> writer(fullFName);

    // Solver interface header and types check
    writer.writeHeader(classHeader(),classVersion());
    writer.writeTypeSizeCheck();    

    // Solver interface parameters header
    writer.writeHeader(params.classHeader(),params.classVersion());
   
    // mesh parameters header
    writer.writeHeader(meshParams.classHeader(),meshParams.classVersion());
    
    // Solver implementation header
    writer.writeHeader(implementation->classHeader(),implementation->classVersion());
  
    //parameters and data
    manager.write(&writer);  
    writeNonParameters(&writer);  

    // Now we can serialize the mesh and implementation
    mesh->write(&writer);    
    implementation->onWrite(&writer); 
    double elapsed = writeRestartTimer->stop();
    glb::console->print<LOG_STD>("done in %.2fs.\n",elapsed);
  }

  void setupNonParameters()
  {
    if (reader != NULL)
      {	
	reader->read(&curStep);
	reader->read(&curTime);
	reader->read(&deltaT);
	reader->read(&nStepsSinceLastCoarsen);
	reader->read(&resimulating);
	restarting=true;
      }
    else
      {
	curStep = 0;
	curTime = params.t0;
	deltaT = params.dt;
	nStepsSinceLastCoarsen=0;
	resimulating=0;
	restarting=false;
      }    
    startStep=curStep;
  }

  template <class W>
  void writeNonParameters(W *writer) const
  {
    writer->write(&curStep);
    writer->write(&curTime);  
    writer->write(&deltaT);
    writer->write(&nStepsSinceLastCoarsen);
    writer->write(&resimulating);
  }

  int checkSignals() const
  {
    //static const unsigned int signal_STOP=(1<<0);
    //static const unsigned int signal_DUMP=(1<<1); 
    int signal=0;

    if (mpiCom->rank()==0)
      {
	if (myIO::exists(params.outputDir+params.stopSignalFileName)) 
	  {
	    myIO::rmFile(params.outputDir+params.stopSignalFileName);
	    signal|=signal_STOP;
	  }
	/*
	if (myIO::exists(params.outputDir+params.dumpSignalFileName)) 
	  {
	    myIO::rmFile(params.outputDir+params.dumpSignalFileName);
	    signal|=signal_DUMP;
	  }
	*/
	if (myIO::exists(params.outputDir+params.dumpRestartSignalFileName)) 
	  {
	    myIO::rmFile(params.outputDir+params.dumpRestartSignalFileName);
	    signal|=signal_DUMP_RESTART;
	  }
      }

    typename TimerPool::Timer timer;
    timer.start();

    if (mpiCom->size()>1)
      glb::console->printFlush<LOG_STD>("Synchronizing MPI processes ...");
	
    mpiCom->template Bcast<int,true>(&signal,0,1);
    
    double elapsed = timer.stop();
    
    if (mpiCom->size()>1)
      glb::console->printFlush<LOG_STD>("done in %.3lgs.\n",elapsed);   
    
    return signal;
  }

  bool checkNeedToStop(int signal)
  {
    static const unsigned int flag_timeLimit=(1<<0);
    static const unsigned int flag_signal=(1<<1);    
    unsigned int flag=0;

    double stopAt = params.timeLimit*params.timeLimitSafety*3600.0;
    bool timerStop = 
      (params.timeLimit>0)&&
      (globalTimer->totalSpent() > (stopAt-stepTimer->lastSpent())); 
  
    if (timerStop) flag|=flag_timeLimit;  
    //if (myIO::exists(params.outputDir+params.stopSignalFileName)) flag|=flag_signal;  
    if (signal&signal_STOP) flag|=flag_signal;
    /*
    typename TimerPool::Timer timer;
    timer.start();

    if (mpiCom->size()>1)
      glb::console->printFlush<LOG_STD>("Synchronizing MPI processes ...");
	
    flag = mpiCom->max<unsigned int,true>(flag);
    
    double elapsed = timer.stop();
    
    if (mpiCom->size()>1)
      glb::console->printFlush<LOG_STD>("done in %.3lgs.\n",elapsed);   
    */

    if (flag&flag_signal)
      {
	glb::console->print<LOG_WARNING>("Stop signal received, stopping.\n");	
	return true;	
      }

     if (flag&flag_timeLimit)
      {
	glb::console->print<LOG_WARNING>("Time limit reached (%.2gs spent), stopping.\n",
					 globalTimer->totalSpent());
	return true;
      }

    return false;
  }
  
  std::string adaptFileName(const std::string &fname, 
			    int mpiRank=-1, 
			    bool restartSuffix=true)
  {
    std::string newFName=fname;
    
    // Add a suffix corresponding to restart step
    if ((restarting)&&(restartSuffix))
      {
	char buf[1024];
	sprintf(buf,"%s_S%6.6ld",newFName.c_str(),startStep);
	newFName=std::string(buf);
      }

    // Check that the file does not already exists
    FILE *f = fopen(newFName.c_str(),"r");
    char buf[1024];
    int index=1;
    if (f!=NULL)
      {
	index=2;
	do {
	  fclose(f);
	  sprintf(buf,"%s_V%d",newFName.c_str(),index++);
	  f = fopen(buf,"r");
	} while (f!=NULL);
	index-=1;
	//newFName=std::string(buf);	    
      }	
    
    if (mpiRank<0) mpiCom->max(index);
    
    if (index>1) 
      {
	sprintf(buf,"%s_V%d",newFName.c_str(),index);
	newFName=std::string(buf);	    
      }
   
    return newFName;
  }

  static std::string buildHeaderString(const std::string &is)
  {
    std::ostringstream oss;
    int inWord=false;
    int count=0;

    for (long i=0;i<is.length();++i)
      {
	if (is[i]==' ')
	  {
	    if (inWord)
	      {
		oss << "("<< count++ <<")"<<" ";	    
		inWord=false;
	      }
	  }
	else 
	  {
	    if ((!inWord)&&(is[i]!='#')) inWord=true;
	    oss << is[i];
	  }
      }
    if (inWord) oss << "("<< count++ <<")"<<" ";

    return oss.str();
  }

  void dumpTimings()
  {
    static std::string timingsFileName=
      params.outputDir+std::string("timings/")+params.timingsFileName;
    static bool initialized=false;
    FILE *timingsFile;
    /*
    if (mpiCom->rank()!=0)  
      {
	return;
      }
    */       
    
    if (!initialized)
      {	
	if (mpiCom->size()>1)
	  {
	    char tmp[255];
	    sprintf(tmp,"_R%4.4d",mpiCom->rank());
	    timingsFileName = timingsFileName + std::string(tmp);
	  }

	//create the directory
	if (mpiCom->rank()==0)
	  myIO::makeDir(params.outputDir+std::string("timings/"));
	mpiCom->barrier();
	timingsFileName = adaptFileName(timingsFileName);

	// Write header
	timingsFile = fopen(timingsFileName.c_str(),"w");
	if (timingsFile==NULL) 
	  {
	    glb::console->print<LOG_ERROR>
	      ("Opening file %s for writing.\n",
	       timingsFileName.c_str());
	    exit(-1);
	  }

	std::ostringstream oss;
	oss << "#stepIndex time nVertices nSimplices nLocalVertices nLocalSimplices nGhostVertices nGhostSimplices";
	implementation->onWriteTimings(curTime,oss,true);
	oss << " " << glb::timerPool->getTimersName();
	fprintf(timingsFile,"%s\n",buildHeaderString(oss.str()).c_str());

	fclose(timingsFile);
	initialized=true; 
      }

    // Write timings
    timingsFile = fopen(timingsFileName.c_str(),"a");
    if (timingsFile==NULL) 
      {
	glb::console->print<LOG_ERROR>
	  ("Opening file %s for appending.\n",
	   timingsFileName.c_str());	
	exit(-1);
      }
    
    char tmp[256];
    std::ostringstream oss;
    implementation->onWriteTimings(curTime,oss,false);

    sprintf(tmp,"%ld %e %ld %ld %ld %ld %ld %ld",
	    curStep,curTime,
	    mesh->getGlobalNCells(0),mesh->getGlobalNCells(Mesh::NDIM),
	    mesh->getNCells(0),mesh->getNCells(Mesh::NDIM),
	    mesh->getNGhostVertices(),mesh->getNGhostSimplices());

    fprintf(timingsFile,"%s%s %s\n",
	    tmp,oss.str().c_str(),glb::timerPool->getTimings().c_str());
    fclose(timingsFile);   
  }

  void dumpStats()
  {
    static std::string statsFileName=
      params.outputDir+params.statisticsFileName;
    static bool initialized=false;
    FILE *statsFile;

    if (mpiCom->rank()!=0)  return;
    
    if (!initialized)
      {	   
	statsFileName = adaptFileName(params.outputDir+params.statisticsFileName,0);

	// Write header
	statsFile = fopen(statsFileName.c_str(),"w");
	if (statsFile==NULL) 
	  {
	    glb::console->print<LOG_ERROR>
	      ("Opening file %s for writing.\n",
	       statsFileName.c_str());
	    exit(-1);
	  }

	std::ostringstream oss;
	oss << "#stepIndex time nVertices nSimplices ";
	implementation->onWriteStatistics(curTime,oss,true);
	fprintf(statsFile,"%s\n",buildHeaderString(oss.str()).c_str());

	fclose(statsFile);
	initialized=true; 
      }

    // Write stats
    statsFile = fopen(statsFileName.c_str(),"a");
    if (statsFile==NULL) 
      {
	glb::console->print<LOG_ERROR>
	  ("Opening file %s for appending.\n",
	   statsFileName.c_str());	
	exit(-1);
      }
    
    char tmp[256];
    sprintf(tmp,"%ld %e %ld %ld",
	    curStep,curTime,
	    mesh->getGlobalNCells(0),mesh->getGlobalNCells(Mesh::NDIM));
    std::ostringstream oss;
    implementation->onWriteStatistics(curTime,oss,false);
    fprintf(statsFile,"%s%s\n",tmp,oss.str().c_str());
    fclose(statsFile);   
  }

  void solve()
  {   
    double elapsed;
    //bool timingsFirstDump=1;
    //const std::string timingsFileName=params.outputDir+params.timingsFileName;
    /*
    if (mpiCom->rank()==0) 
      {
	timingsFile = fopen(timingsFileName.c_str(),"w");
	if (timingsFile==NULL) 
	  {
	    glb::console->print<LOG_ERROR>
	      ("Opening file %s for writing.\n",
	       (params.outputDir+params.timingsFileName).c_str());
	    exit(-1);
	  }
      }
    */
    bool restarting = (reader!=NULL);
    globalTimer->start();

    buildTimer->start();
    if (restarting) mesh->read(meshParams,reader);
    else implementation->onBuildMesh(mesh,meshParams);
    elapsed = buildTimer->check();
    glb::console->print<LOG_INFO>("Mesh was built in %lg seconds.\n",elapsed);
    
    glb::console->print<LOG_INFO>
      ("Calling '%s' initializer.\n",SolverImplementation::classHeader().c_str());
    glb::console->indent();
    if (restarting) implementation->onRead(reader);
    implementation->onInitialize(mesh,params,curTime,restarting);
    elapsed = buildTimer->check()-elapsed;
    glb::console->unIndent();
    glb::console->print<LOG_INFO>("Done in %lg seconds.\n",elapsed);
    
    glb::console->unIndent();  
    elapsed = buildTimer->stop();
    glb::console->print<LOG_STD>("Solver initialization completed in %lgs.\n",elapsed);
    glb::memoryInspector->report<LOG_STD>();

    if (reader!=NULL) 
      {
	delete reader;
	reader = BinaryReader::nullReader();
      }
    //else writeRestartFile("initialRestart");    
    
    bool mustStop=false;
    long nStepsSinceLastRestart=0; 
    double el[5];
    int signals = checkSignals();

    //glb::console->print<LOG_STD>("Starting solver @t=%lg\n",curTime);
  
    PPROF_START_DEFAULT;

  SolverMainLoop:    
    for (; 
	 ((curTime<params.tEnd+(deltaT/2))&&(!mustStop)); 
	 curTime+=deltaT, curStep++)
      {		

	if (((nStepsSinceLastRestart>=params.restartEvery)&&(params.restartEvery>0))||
	    (signals & signal_DUMP_RESTART))
	  {
	    bool fromSignal=signals&signal_DUMP_RESTART;
	    char fname[256];
	    sprintf(fname,"step_%6.6ld",curStep);
	    //mpiCom->waitMyTurn(1);	    
	    writeRestartFile(fname, fromSignal);	    
	    //mpiCom->finishedMyTurn(1);
	    
	    if ((nStepsSinceLastRestart>=params.restartEvery)&&(params.restartEvery>0))
	      nStepsSinceLastRestart=1;
	    else 
	      nStepsSinceLastRestart++;
	  }
	else nStepsSinceLastRestart++;
	
	glb::console->print<LOG_STD>("Starting time step %ld (t=%lg).\n",curStep,curTime);
	glb::console->indent();

	stepTimer->start();		
	// hlp::MakeConst is to ensure that the implementation will not change the value
	updateTimer->start();
	implementation->onNewTimeStep
	  (hlp::makeConst(curStep),hlp::makeConst(curTime),deltaT);
	el[0]=updateTimer->stop();

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_START("Advance");}

	advanceTimer->start();
	implementation->onAdvance();
	el[1]=advanceTimer->stop();

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_STOP("Advance");}

	if (glb::console->willPrint<LOG_INFO>())
	  {
	    glb::console->print<LOG_INFO>("Updating unstructured mesh:\n");
	    glb::console->indent();
	  }

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_START("Coarsen");}

	coarsenTimer->start();
	long nCoarsened=coarsen();
	implementation->afterCoarsen(nCoarsened);
	el[2]=coarsenTimer->stop();
	
#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_STOP("Coarsen");}

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_START("Refine");}

	refineTimer->start();
	long nRefined=refine();
	implementation->afterRefine(nRefined);
	el[3]=refineTimer->stop();

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_STOP("Refine");}

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_START("Repart");}

	repartTimer->start();
	bool repartStatus = repart();
	implementation->afterRepart(repartStatus);
	el[4]=repartTimer->stop();

#pragma omp parallel num_threads(glb::num_omp_threads)
	{LIKWID_MARKER_STOP("Repart");}	

	elapsed = stepTimer->stop();

	dumpTimings();
	dumpStats();

	if (glb::console->willPrint<LOG_INFO>())
	  {
	    glb::console->unIndent();
	    glb::console->print<LOG_INFO>("Done in %lgs\n",el[2]+el[3]+el[4]);
	  }	

	glb::console->unIndent();	
	
	if (glb::console->willPrint<LOG_INFO>())
	  {
	    // glb::console->print<LOG_INFO>("All done in %lgs (%.3lg/%.3lg/%.3lg/%.3lg/%.3lg).\n",
	    // 				  elapsed,el[0],el[1],el[2],el[3],el[4]);
	    glb::console->print<LOG_INFO>("All done in %lgs (%.2f / %.2f / %.2f / %.2f).\n",
					  elapsed,el[1],el[2],el[3],el[4]);
	  }
	else
	  glb::console->print<LOG_STD>("All done in %lgs.\n",elapsed);
	
	if (glb::console->willPrint<LOG_PEDANTIC>())
	  {
	    // NB: This requires a global synchronization, so it is probably a big waste
	    // of time !
	    // Should print the local values for each process instead, and forget about the
	    // imbalance factor
	    // glb::console->print<LOG_PEDANTIC>("The mesh now has %ld vertices and %ld simplices\n",
	    // 				      mesh->getGlobalNCells(0),
	    // 				      mesh->getGlobalNCells(M::NDIM));
	    glb::console->
	      print<LOG_PEDANTIC>("The local mesh now has %ld(+%ldg/s) vertices and %ld(+%ldg/s) simplices\n",
				  mesh->getNCells(0),
				  mesh->getNCellsTotal(0)-mesh->getNCells(0),
				  mesh->getNCells(Mesh::NDIM),
				  mesh->getNCellsTotal(Mesh::NDIM)-mesh->getNCells(Mesh::NDIM));

	    glb::memoryInspector->report<LOG_PEDANTIC>();

	    // glb::console->print<LOG_PEDANTIC>("New imbalance factor is %.3f.\n",
	    // 				      mesh->getLoadImbalanceFactor());
	  }

	signals = checkSignals();
	// Synchronize processes and check whether we need to stop
	// do not use 'break' here so that counters are correctly incremented, 
	mustStop = checkNeedToStop(signals);
	/*
	if (mpiCom->size()>1)
	  glb::console->printFlush<LOG_STD>("Synchronizing MPI processes ...");	
	int stopSignal = mpiCom->max<int,true>((int)myIO::exists());
	if (mpiCom->size()>1)
	  glb::console->printFlush<LOG_STD>("done in %.3lgs.\n",);
	
	
	if (stopSignal)
	  {
	    glb::console->print<LOG_STD>("Stop signal received, stopping.\n");
	    // do not use 'break' here so that counters are correctly incremented, 
	    mustStop=true;
	  }

	if (needToStop()) 
	  {
	    glb::console->print<LOG_STD>("Time limit reached (%.2gs spent), stopping.\n",
					 globalTimer->totalSpent());
	    // do not use 'break' here so that counters are correctly incremented, 
	    mustStop=true;
	  }
	*/

      } // SolverMainLoop
     
    
    if (true)
      {
	char fname[256];
	sprintf(fname,"step_%6.6ld",curStep);
	//mpiCom->waitMyTurn(1);
	writeRestartFile(fname);
	//mpiCom->finishedMyTurn(1);	
      }
 
    mpiCom->barrier();
    glb::memoryInspector->report<LOG_INFO>();
    elapsed=stepTimer->totalSpent();
    glb::console->print<LOG_STD>("Finished computing in %lgs.\n",elapsed);

    PPROF_STOP;
    LIKWID_MARKER_CLOSE;

    if (glb::console->willPrint<LOG_INFO>())
      {
	el[0] = updateTimer->totalSpent();
	el[1] = advanceTimer->totalSpent();
	el[2] = coarsenTimer->totalSpent();
	el[3] = refineTimer->totalSpent();
	el[4] = repartTimer->totalSpent();
	
	
	glb::console->indent();
	// glb::console->print<LOG_INFO>("=> Spent %lgs updating (%.1f%%)\n",
	// 			 el[0],el[0]/elapsed*100);
	glb::console->print<LOG_INFO>("=> Spent %lgs advancing (%.1f%%)\n",
				      el[1],el[1]/elapsed*100);
	glb::console->print<LOG_INFO>("=> Spent %lgs coarsening (%.1f%%)\n",
				      el[2],el[2]/elapsed*100);
	glb::console->print<LOG_INFO>("=> Spent %lgs refining (%.1f%%)\n",
				      el[3],el[3]/elapsed*100);	
	glb::console->print<LOG_INFO>("=> Spent %lgs repartitionning (%.1f%%)\n",
				      el[4],el[4]/elapsed*100);
	
	double global[3],local[3],barrier[3],total[3];
	int rk[2];
	const char *text[3]={"Minimum","Maximum","Average"};
	mpiCom->reportTimeSpent(global,local,barrier,total,rk);
	if (total[0]>0)
	  {
	    int i;
	    for (i=0;i<2;++i)
	      {
		glb::console->print<LOG_INFO>("=> %s time spent communicating is %lgs (%.1f%, rk=%d) : [%.1f%% global, %.1f%% local, %.1f%% barriers]\n", 
					      text[i],total[i],total[i]/elapsed*100,rk[i],
					      global[i]/total[i]*100,local[i]/total[i]*100,
					      barrier[i]/total[i]*100);
	      }
	    glb::console->print<LOG_INFO>("=> %s time spent communicating is %lgs (%.1f%%) : [%.1f%% global, %.1f%% local, %.1f%% barriers]\n", 
					  text[i],total[i],total[i]/elapsed*100,
					  global[i]/total[i]*100,local[i]/total[i]*100,
					  barrier[i]/total[i]*100);
	  }
	glb::console->unIndent();
      }

    if (params.resimulate && (!resimulating) && (!mustStop))
      {
	glb::console->print<LOG_STD>("\n");
	glb::console->print<LOG_STD>("First pass finished, starting resimulation.\n");
	glb::console->print<LOG_STD>("\n");
	
	setupNonParameters();
	resimulating = true;
	implementation->onResimulate();
	goto SolverMainLoop;
      }
    
    //if (mpiCom->rank()==0) fclose(timingsFile);   
  }
 
  long refine(hlp::IsEnabled)
  {    
    if (params.noRefine)
      return 0;
    else
      {
	//glb::console->printFlush<LOG_STD>("Refining took  ");  
	return mesh->refine(implementation);    	
      }
  }

  long refine(hlp::IsDisabled)
  {
    return 0;
  }

  long refine()
  {
    return refine(RefinementStatus());
  }

  long coarsen(hlp::IsEnabled)
  {    
    long result=0;
    if (!params.noCoarsen)
      {
	if (nStepsSinceLastCoarsen == params.coarsenEvery)
	  {
	    nStepsSinceLastCoarsen=1;
	    result=mesh->coarsen(implementation);    
	  }
	else nStepsSinceLastCoarsen++;	
      }
    return result;
  }

  long coarsen(hlp::IsDisabled)
  {
    return 0;
  }

  long coarsen()
  {
    return coarsen(CoarseningStatus());
  }

  bool repart()
  {
    //double weight=-1; // Use default weight
    bool force=false;
    double weight=implementation->onCheckRepartLocalWeight(stepTimer->check(), force);
    return mesh->repart(weight,force);
  }
 
private:
  Mesh *mesh;
  SolverImplementation *implementation;  
  BinaryReader *reader;
  ParamsManager manager;
 
  MpiCommunication *mpiCom;

  Params params;
  MeshParams meshParams;

  bool fromRestartFile;
  
  long curStep;
  double curTime;
  double deltaT;  
  long nStepsSinceLastCoarsen;
  long resimulating;

  bool restarting;
  long startStep;

  float version;

  typename TimerPool::Timer *globalTimer;
  typename TimerPool::Timer *buildTimer;
  typename TimerPool::Timer *stepTimer;
  typename TimerPool::Timer *advanceTimer;
  typename TimerPool::Timer *coarsenTimer;
  typename TimerPool::Timer *refineTimer;
  typename TimerPool::Timer *repartTimer;  
  typename TimerPool::Timer *updateTimer;  
  typename TimerPool::Timer *writeRestartTimer;  
};

} // namespace slv

#include "../internal/namespace.footer"
#endif
