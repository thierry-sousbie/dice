#ifndef __COLDICE_CONFIG_H__
#define __COLDICE_CONFIG_H__

#include <dice/tools/wrappers/boostMultiprecisionFloat128.hxx>
#include <dice/tools/helpers/helpers_macros.hxx>

// Because some compilers do not accept spaces in -D[...] macro definition options ...
typedef long double long_double;
typedef long double longdouble;

// Floating point types used for projection
#ifndef D_PROJECTION_FLOAT_TYPE
#define D_PROJECTION_FLOAT_TYPE dice::DoubleDouble
#endif

// High resolution floating point types used for projection
#ifndef D_PROJECTION_HR_FLOAT_TYPE
#define D_PROJECTION_HR_FLOAT_TYPE dice::DoubleDouble
#endif

// Default Number of dimensions (in real space)
#ifndef D_DIMS_COUNT
#define D_DIMS_COUNT 2
#endif

// Root nodes level for AMR grid
#ifndef D_AMR_ROOT_LEVEL
#define D_AMR_ROOT_LEVEL 6
#endif

// boundary conditions (can be PERIODIC or BOXED)
#ifndef D_BOUNDARY_TYPE
#define D_BOUNDARY_TYPE PERIODIC
#endif

// Maximum number of threads for drift (8 should be already more than enough)
#ifndef D_N_THREADS_MAX_DRIFT
#define D_N_THREADS_MAX_DRIFT 8
#endif

// Enable/Disable accuracy checking for projection
#ifndef D_ENABLE_ACCURACY_CHECKING
#define D_ENABLE_ACCURACY_CHECKING false
#endif

// Define if you don't need simplices tracers (you most probably don't)
#define NO_SIMPLEX_TRACERS

// define to allow using more memory to improve the speed of some routines
//#define D_TRADE_MEMORY_FOR_SPEED

// define to differenciate potential via FFT instead of using finite differences
//#define D_FFT_POISSON_SOLVER_COMPUTE_DISPLACEMENT

// define to ignore tracers when refining the mesh 
// => split the longest edge and refine using midPoints but tracers are still advected
//#define D_DONT_USE_TRACERS

namespace cfg {
  
  template <class LOG,class C>
    void printDefines(C *console)
    {

      console->template print<LOG>
	("Machine configuration:\n");
      console->indent();
      console->template print<LOG>
	("sizeof(long) is %ld bytes.\n",sizeof(long));  
      console->template print<LOG>
	("sizeof(long double) is %ld bytes (%d significant bits).\n",
	 sizeof(long double),std::numeric_limits<long double>::digits); 
      if (dice::Float128OrMore_IsQuadruplePrecision::value)
	console->template print<LOG>
	  ("Quadruple precision 'float128' type is available as dice::Float128OrMore (%d significant bits).\n",
	   std::numeric_limits<dice::Float128OrMore>::digits); 
      else
	{
	  console->template print<LOG>
	    ("Quadruple precision 'float128' type is NOT available: using gmp_float as dice::Float128OrMore instead (%d significant bits).\n",std::numeric_limits<dice::Float128OrMore>::digits);
	}

      console->unIndent();
      console->template print<LOG>("\n");
      console->template print<LOG>
	("Solver was compiled with the following configuration:\n");
      console->indent();
      console->template print<LOG>
	("DIMS_COUNT = %d\n",D_DIMS_COUNT);
      console->template print<LOG>
	("AMR_ROOT_LEVEL = %d\n",D_AMR_ROOT_LEVEL);
      console->template print<LOG>
	("N_THREADS_MAX_DRIFT = %d\n",D_N_THREADS_MAX_DRIFT);
      console->template print<LOG>
	("ENABLE_ACCURACY_CHECKING = %d\n",D_ENABLE_ACCURACY_CHECKING);
      console->template print<LOG>
	("PROJECTION_FLOAT_TYPE = %s\n", STRINGIFY(D_PROJECTION_FLOAT_TYPE) );
      console->template print<LOG>
	 ("PROJECTION_HR_FLOAT_TYPE = %s\n", STRINGIFY(D_PROJECTION_HR_FLOAT_TYPE) );
      /* console->template print<LOG> */
      /* 	("BOUNDARY_TYPE = %s\n",STRINGIFY(D_BOUNDARY_TYPE)); */
      console->unIndent();
      console->template print<LOG>("\n");
    }
}


#endif
