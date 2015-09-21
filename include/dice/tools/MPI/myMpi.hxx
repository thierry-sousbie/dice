#ifndef __MY_MPI_HXX__
#define __MY_MPI_HXX__

#ifdef USE_MPI
#include <mpi.h>

#ifndef MPID_TAG_UB
#define MPID_TAG_UB (1<<((sizeof(int)*8) -2)) // maximum tag value
#endif

#include "../../internal/namespace.header"

template <typename T>
struct MPI_Type
{
  static inline MPI_Datatype get(){return MPI_DATATYPE_NULL;}
};

#define ASSOCIATE_MPI_TYPE(T,M)			\
  template <>					\
  struct MPI_Type<T> {				\
    static inline MPI_Datatype get(){return M;}	\
  };						\

 
ASSOCIATE_MPI_TYPE(char              ,MPI_CHAR)
ASSOCIATE_MPI_TYPE(unsigned char     ,MPI_UNSIGNED_CHAR)
ASSOCIATE_MPI_TYPE(short             ,MPI_SHORT)
ASSOCIATE_MPI_TYPE(unsigned short    ,MPI_UNSIGNED_SHORT)
ASSOCIATE_MPI_TYPE(int               ,MPI_INT)
ASSOCIATE_MPI_TYPE(unsigned int      ,MPI_UNSIGNED)
ASSOCIATE_MPI_TYPE(long              ,MPI_LONG)
ASSOCIATE_MPI_TYPE(unsigned long     ,MPI_UNSIGNED_LONG)
ASSOCIATE_MPI_TYPE(long long         ,MPI_LONG_LONG)
ASSOCIATE_MPI_TYPE(unsigned long long,MPI_UNSIGNED_LONG_LONG)
ASSOCIATE_MPI_TYPE(float             ,MPI_FLOAT)
ASSOCIATE_MPI_TYPE(double            ,MPI_DOUBLE)
ASSOCIATE_MPI_TYPE(long double       ,MPI_LONG_DOUBLE)

/*
ASSOCIATE_MPI_TYPE(int8_t  ,MPI_INT8_T)
ASSOCIATE_MPI_TYPE(int16_t ,MPI_INT16_T)
ASSOCIATE_MPI_TYPE(int32_t ,MPI_INT32_T)
ASSOCIATE_MPI_TYPE(int64_t ,MPI_INT64_T)

ASSOCIATE_MPI_TYPE(uint8_t  ,MPI_UINT8_T)
ASSOCIATE_MPI_TYPE(uint16_t ,MPI_UINT16_T)
ASSOCIATE_MPI_TYPE(uint32_t ,MPI_UINT32_T)
ASSOCIATE_MPI_TYPE(uint64_t ,MPI_UINT64_T)
*/
#include "../../internal/namespace.footer"
#else // we don't have MPI

#ifndef MPID_TAG_UB
#define MPID_TAG_UB (1<<((sizeof(int)*8) -2))
#endif
#define MPI_ANY_TAG 0

#ifndef MPI_ORDER_FORTRAN
#define MPI_ORDER_FORTRAN 1
#endif
#ifndef MPI_ORDER_C
#define MPI_ORDER_C 2
#endif

typedef void* MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Request;
typedef void* MPI_Op;
typedef int MPI_Aint;
struct MPI_Status {static const int MPI_SOURCE=0;};

MPI_Op MPI_MAX = 0;
MPI_Op MPI_MIN = 0;
MPI_Op MPI_SUM = 0;
MPI_Datatype MPI_BYTE = 1;
MPI_Datatype MPI_DATATYPE_NULL = 0;

#include "../../internal/namespace.header"

template <typename T>
struct MPI_Type
{
  static inline MPI_Datatype get(){return MPI_DATATYPE_NULL;}
};

#include "../../internal/namespace.footer"
  
int MPI_Type_free(MPI_Datatype *datatype) {return 0;}
void MPI_Type_create_struct(long,int*,MPI_Aint*,MPI_Datatype*,MPI_Datatype*) {}
void MPI_Type_commit(MPI_Datatype*) {}
int MPI_Type_create_subarray(int ndims,
			     int array_of_sizes[],
			     int array_of_subsizes[],
			     int array_of_starts[],
			     int order,
			     MPI_Datatype oldtype,
			     MPI_Datatype *newtype) {return 0;}
int MPI_Type_contiguous(int count,
			MPI_Datatype old_type,
			MPI_Datatype *new_type_p) {return 0;}

/*
// This MUST be called once from a single thread first for initialization
double MPI_Wtime()
{
static time_t refTime;
  static bool init=true;

  if (init)
    {
      init=false;
      time(&refTime);    
    }

  time_t now;
  time(&now);
  return difftime(now,refTime);
}
*/
#ifdef __MACH__ // This is OSX
#include <mach/mach.h>
#include <mach/mach_time.h>

double MPI_Wtime()
{   
  static double timeConvert = 0.0;
  static double ref;
  if ( timeConvert == 0.0 )
    {
      mach_timebase_info_data_t timeBase;
      (void)mach_timebase_info( &timeBase );
      timeConvert = (double)timeBase.numer /
	(double)timeBase.denom /
	1000000000.0;
      ref=(double)mach_absolute_time( )*timeConvert;
    }
  return (double)mach_absolute_time( )*timeConvert - ref;
}

#else // NOT OSX

#include <time.h>

// This MUST be called once from a single thread first for initialization
double MPI_Wtime()
{ 
  static timespec ref;
  static bool init=true;

  if (init)
    {
      init=false;      
      clock_gettime(CLOCK_REALTIME, &ref);     
    }

  timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  return double(now.tv_sec-ref.tv_sec) + 1.E-9L * double(now.tv_nsec-ref.tv_nsec);
}

#endif // OSX ?

#endif // MPI ?

#endif // HEADER
