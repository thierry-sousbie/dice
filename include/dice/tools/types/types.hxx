#ifndef __MY_TYPES_HXX__
#define __MY_TYPES_HXX__

#include "../../dice_globals.hxx"

#include "../../tools/helpers/helpers.hxx"


#include "../../internal/namespace.header"

#ifdef USELONGINT
typedef long int INT;
typedef unsigned long int UINT;
#else
typedef int INT;
typedef unsigned int UINT;
#endif

#ifdef USESIMPLEPRECISION
typedef float FLOAT;
#else
typedef double FLOAT;
#endif

typedef long long int Long64;
typedef unsigned long long int ULong64;
/*
struct CheckTypeSize
{
  template <class T, int size>
  static int warnTypeSize(const char *name)
  {
    if (!hlp::IsEqual<sizeof(T),size>::result)
      {
	if (sizeof(T)<size)
	  {
	    PRINT_SRC_INFO(LOG_WARNING);
	    console->print<LOG_WARNING>("Type %s is smaller than expected. (%d bytes, should be %d)\n",
				      name,sizeof(T),size);
	    console->print<LOG_WARNING>("You should edit file 'types.hxx'.\n");	    
	    return -1;
	  }
	else
	  {
	    PRINT_SRC_INFO(LOG_WARNING);
	    console->print<LOG_WARNING>("Type %s is bigger than expected. (%d bytes, should be %d)\n",
					name,sizeof(T),size);
	    console->print<LOG_WARNING>("You should edit file 'types.hxx'.\n");
	    return +1;
	  }
      }
  }

};
*/

#include "../../internal/namespace.footer"
#endif
