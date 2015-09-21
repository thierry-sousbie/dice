#ifndef __LOG_TRAITS_HXX__
#define __LOG_TRAITS_HXX__

#include <stdio.h>


#include "../../internal/namespace.header"

struct LOG_STD
{
  static const int verboseLevel = 1;
  static const bool mpiMainRankOnly =true;
  static const int nOutput = 1;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_STD::output[]={stdout};
const char LOG_STD::header[]="";

struct LOG_STD_ALL : LOG_STD
{
  static const bool mpiMainRankOnly =false;
};

struct LOG_INFO
{
  static const int verboseLevel = 2;
  static const bool mpiMainRankOnly =true;
  static const int nOutput = 1;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_INFO::output[]={stdout};
const char LOG_INFO::header[]="";

struct LOG_INFO_ALL : LOG_INFO
{
  static const bool mpiMainRankOnly =false;
};

struct LOG_PEDANTIC
{
  static const int verboseLevel = 3;
  static const bool mpiMainRankOnly =true;
  static const int nOutput = 1;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_PEDANTIC::output[]={stdout};
const char LOG_PEDANTIC::header[]="";

struct LOG_PEDANTIC_ALL : LOG_PEDANTIC
{
  static const bool mpiMainRankOnly =false;
};

struct LOG_DEBUG
{
  static const int verboseLevel = 4;
  static const bool mpiMainRankOnly =false;
  static const int nOutput = 1;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_DEBUG::output[]={stdout};
const char LOG_DEBUG::header[]="[DEBUG] ";

struct LOG_WARNING
{
  static const int verboseLevel = 0;
  static const bool mpiMainRankOnly =false;
  static const int nOutput = 2;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_WARNING::output[]={stdout,stderr};
const char LOG_WARNING::header[]="[WARNING] ";

struct LOG_WARNING_SINGLE
{
  static const int verboseLevel = 0;
  static const bool mpiMainRankOnly =true;
  static const int nOutput = 2;
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_WARNING_SINGLE::output[]={stdout,stderr};
const char LOG_WARNING_SINGLE::header[]="[WARNING] ";

struct LOG_ERROR
{
  static const int verboseLevel = 0;
  static const bool mpiMainRankOnly =false;
  static const int nOutput = 2; 
  static const char header[];
  static FILE* output[nOutput];
};
FILE* LOG_ERROR::output[]={stdout,stderr};
const char LOG_ERROR::header[]="<<ERROR>> ";

#include "../../internal/namespace.footer"
#endif
