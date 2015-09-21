#include <string>
#include <stdlib.h>
#include <string.h>
#include <sstream> 

#include "dice-config.hxx"

void usage(char *exe)
{
  //printf("\n");
  printf("Usage: %s [-cxxflags] [-libs] [-defines] [-includes] [-all] [-all_noflags]\n",exe);
}

int main(int argc, char **argv)
{
  if (argc==1) usage(argv[0]);
  //else printf("\n");
  int nFound=0;
  int errors=0;
  for (int i=1;i<argc;)
    {
      std::string arg=argv[i++];
      int oldFound=nFound;
      int all=arg==std::string("-all");
      int allButFlags=arg==std::string("-all_noflags");
      if (allButFlags) all=1;
      
      if ((arg==std::string("-compiler"))||all) {printf(" %s ",COMPILER);nFound++;}
      if (((arg==std::string("-flags"))||all)&&(!allButFlags)) {printf(" %s ",CXX_FLAGS);nFound++;}
      if ((arg==std::string("-includes"))||all) {printf(" %s ",INCLUDES);nFound++;}
      if ((arg==std::string("-defines"))||all) {printf(" %s ",DEFINES);nFound++;}
      if ((arg==std::string("-libs"))||all) {printf(" %s ",LIBRARIES);nFound++;}
      
      if (oldFound==nFound) 
	{
	  printf("\nERROR: unknown parameter: '%s'\n",argv[i-1]);
	  errors=1;
	}
    }  
  printf("\n");
  
  if (errors) 
    {
      usage(argv[0]);
      return -1;
    }

  return 0;
}
