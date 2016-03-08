#ifdef GLOBAL_DEFINITION
#define GLOBAL
#else 
#define GLOBAL extern
#endif


#include "./internal/namespace.header"

class MpiCommunication;
class Console;
class ParamsParser;
class MemoryInspector;
class TimerPool;

namespace glb {
  GLOBAL MpiCommunication *mpiComWorld;
  GLOBAL Console *console;
  GLOBAL ParamsParser *pParser;
  GLOBAL ParamsParser *dummyPParser;
  GLOBAL MemoryInspector *memoryInspector;
  GLOBAL TimerPool *timerPool;
  GLOBAL int debug;
  GLOBAL int num_omp_threads;
  GLOBAL int global_stop_signal;

  GLOBAL int debugMeshRefine;
}

#include "./internal/namespace.footer"

#undef GLOBAL
#undef GLOBAL_DEFINITION

#ifndef __GLOBAL_DEFINITION_HXX__
#define __GLOBAL_DEFINITION_HXX__

#define PRINT_SRC_INFO(logtype) glb::console->printSrcInfo<logtype>(__PRETTY_FUNCTION__,__FUNCTION__,__FILE__,__LINE__);

#include "./tools/MPI/myMpi.hxx" // Need to include before timerPool with intelMPI ...
#include "./tools/time/timerPool.hxx"
#include "./tools/MPI/mpiCommunication.hxx"
#include "./tools/IO/console.hxx"
#include "./tools/IO/paramsParser.hxx"
#include "./tools/OMP/openMP_interface.hxx"
#include "./tools/memory/memoryInspector.hxx"

#ifdef HAVE_EIGEN3
#ifdef NDEBUG
#define EIGEN_NO_DEBUG
#endif
#endif

#ifdef USE_OPENMP
#ifdef HAVE_EIGEN3
// gcc complains about deprecated declarations in eigen ...
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// eigen needs to be initialized when multithreading is enabled !
#include <Eigen/Core>
#pragma GCC diagnostic pop
#endif
#endif

#include "./internal/namespace.header"

/** \brief sets the default number of threads for the entire library
 */
void setGlobalNumThreads(int nThreads)
{
  omp_set_num_threads(nThreads);
  glb::num_omp_threads=nThreads;
  glb::console->print<LOG_INFO>("Changing default number of openMP threads to %d.\n",
				glb::num_omp_threads);
}

struct GlobalObjects
{ 
  static void init(int *argc, char ***argv,
		   bool enableParamsParser=true,
		   bool enableGlobalLog=true,
		   bool initializeMpi=true)
  {
    glb::num_omp_threads=1;
    glb::debug=0;  
    // PAY ATTENTION TO THE ORDER HERE !!!!
    glb::timerPool = new TimerPool();

    if (initializeMpi)
      glb::mpiComWorld = new MpiCommunication(argc,argv);
    else
      glb::mpiComWorld = new MpiCommunication();

    glb::console = new Console(3); // parameter is the default verbose level
    glb::memoryInspector = new MemoryInspector();
    if (enableParamsParser)
      glb::pParser = new ParamsParser(*argc,*argv);
    else
      glb::pParser = new ParamsParser();
    glb::dummyPParser = new ParamsParser;  
    glb::global_stop_signal=0;

    glb::console->setVerboseLevel(glb::pParser->
				  get<>("verbose",
					ParamsParser::defaultCategory(),
					glb::console->getVerboseLevel(),
					"Verbosity level (from 0=quiet to 4=debug)")
				  );
    
    int noGlobalLog=glb::pParser->get<>("noGlobalLog",
					ParamsParser::defaultCategory(),
					enableGlobalLog?0:1,
					"Set to prevent creation of a global log");
    std::string globalLogDir("");
    globalLogDir=glb::pParser->get<>("globalLogDir",
				     ParamsParser::defaultCategory(),
				     globalLogDir,
				     "The directory where the global log is created");
    if (globalLogDir.length()>0)
      if (globalLogDir[globalLogDir.length()-1] != '/')
	globalLogDir+="/";
    
    std::string globalLogFName("globalLog");
    globalLogFName=glb::pParser->get<>("globalLogFName",
				       ParamsParser::defaultCategory(),
				       globalLogFName,
				       "The name of the global log file");
    
    if (!noGlobalLog) 
      {
	char fname[1024];
	sprintf(fname,"%s%s_%6.6d.log",
		globalLogDir.c_str(),
		globalLogFName.c_str(),
		glb::mpiComWorld->rank());
	glb::console->outputToFile(fname);
      }

    glb::debug=glb::pParser->get<>("debug",
				   ParamsParser::defaultCategory(),
				   glb::debug,
				   "Set for debug mode (value 0 to 2)");
      
    glb::num_omp_threads = glb::pParser->get<>("threads_per_node",
					       ParamsParser::defaultCategory(),
					       -1,
					       "The number of openMP threads to use per MPI node (auto set if <0)");
    bool prtThreads=true;
    if (glb::num_omp_threads<0)
      {
	glb::num_omp_threads = omp_get_max_threads();
	/*
#pragma omp parallel   
	if (omp_get_thread_num()==0) 
	  glb::num_omp_threads=omp_get_num_threads();	
	*/
      }     
    omp_set_num_threads(glb::num_omp_threads);

    int actualNumThreads=omp_get_max_threads();
    /*
#pragma omp parallel   
    if (omp_get_thread_num()==0) 
      actualNumThreads=omp_get_num_threads();
    */    

    if (actualNumThreads!=glb::num_omp_threads)
      {
	glb::console->print<LOG_WARNING>("Could not allocate the requested number of openMP threads (=%d).\n",glb::num_omp_threads);
#ifndef USE_OPENMP
	glb::console->print<LOG_WARNING>("To use openMP, enable it in cmake (-DUSE_OPENMP=true) and recompile.\n");
#endif
	glb::num_omp_threads = actualNumThreads;
	glb::console->print<LOG_WARNING>("Default number of openMP threads set to %d.\n",glb::num_omp_threads);
	
	prtThreads=false;
      }

    if ((prtThreads)&&(glb::num_omp_threads>=1))
      glb::console->print<LOG_INFO>("Default number of openMP threads set to %d.\n",
				    glb::num_omp_threads);
#ifdef HAVE_EIGEN3
#ifdef USE_OPENMP
    // Initialize EIGEN if necessary
    Eigen::initParallel();
    Eigen::setNbThreads(1);
#endif
#endif
    
    glb::debugMeshRefine=0;
    glb::debugMeshRefine=glb::pParser->get<>("meshRefine",
					     "glbDebug",
					     glb::debugMeshRefine,
					     "Set to debug mesh refinement.");
    
  }

  static void init()
  {
    int argc=0;
    char **argv=NULL;
    init(&argc,&argv,false);
  }

  static void cleanUp()
  {    
    delete glb::console;
    delete glb::pParser;
    delete glb::dummyPParser;
    delete glb::memoryInspector;

    glb::mpiComWorld->barrier();

    delete glb::timerPool;
    delete glb::mpiComWorld;
  }  
};

#include "./internal/namespace.footer"
#endif
