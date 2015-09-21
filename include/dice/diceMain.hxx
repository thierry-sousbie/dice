#ifndef __DICE_MAIN_HXX__
#define __DICE_MAIN_HXX__

// This should be included in the cpp file where main() function is defined.

#define GLOBAL_DEFINITION
#include "dice_globals.hxx"
#undef GLOBAL_DEFINITION

#include "./internal/namespace.header"

void initialize(int *argc, char ***argv, bool initializeMpi=true)
{
  GlobalObjects::init(argc,argv,initializeMpi);
}

void finalize()
{  
  GlobalObjects::cleanUp();
}

#include "./internal/namespace.footer"
#endif
