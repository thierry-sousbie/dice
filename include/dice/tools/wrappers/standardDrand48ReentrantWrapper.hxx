#ifndef __STANDARD_DRAND48_REENTRANT_WRAPPER_HXX__
#define __STANDARD_DRAND48_REENTRANT_WRAPPER_HXX__

#include <cstdlib>
#include "../helpers/helpers.hxx"

/**
 * @file 
 * @brief A wrapper to drand48_r and srand48_r on systems that have it and a replacement
 * using erand48 on systems that do not.
 * @author Thierry Sousbie 
 */

// That's ugly, but drand48_r is in :: and I do not know how to import it to another namespace
struct drand48_r_does_not_exist{};
drand48_r_does_not_exist drand48_r(...);
struct drand48_data;

#include "../../internal/namespace.header"

namespace internal {
  
  // Check if the free function int drand48_r(drand48_data *, double *) exists
  // FIXME: I DO NOT KNOW HOW TO IMPORT THE GLOBAL NAMESPACE !!!!!
  // -> this does not work
  /*
  namespace stupid {
    using namespace std; // should import :: here 
    struct does_not_exist{};
    does_not_exist drand48_r(...);
  }
  
  drand48_data dummy_drand48_data;
  double dummy_double;

  typedef hlp::SameType<decltype(drand48_r(&dummy_drand48_data,&dummy_double)),
			drand48_r_does_not_exist> no_drand48_r; 
*/
  typedef hlp::HasDestructor<drand48_data> have_drand48_r;

}

// no drand48_r !
template <bool no_drand48>
struct DRand48_rWapperT {
  static const int FROM_DRAND48_R = 0;
  struct RandData {unsigned short xsubi[3];};
  
  static int drand48_r(RandData *buffer, double *result)
  {    
    *result=::erand48(buffer->xsubi);
    return 0;
  }
  
  static int srand48_r(long int seedval, RandData *buffer)
  {
    unsigned char *b=reinterpret_cast<unsigned char*>(buffer);
    for (int i=0;i<sizeof(RandData);++i)
      {
	*b=seedval&(0xFF);
	seedval>>=8;
	b++;
      }
    return 0;
  }
};

//with drand_48r
/*
template <>
struct DRand48_rWapperT<false>{
  static const int FROM_DRAND48_R = 1;  
  typedef drand48_data RandData;

  static int drand48_r(RandData *buffer, double *result)
  {
    return ::drand48_r(buffer,result);
  }

  static int srand48_r(long int seedval,RandData *buffer)
  {
    return ::srand48_r(seedval,buffer);
  }
};
*/

// wraps to drand48_r ans srand48_r if they exist, or to erand48 otherwise
//typedef DRand48_rWapperT< internal::no_drand48_r::value >  DRand48_rWapper;
typedef DRand48_rWapperT< !internal::have_drand48_r::value >  DRand48_rWapper;

#include "../../internal/namespace.footer"
#endif
