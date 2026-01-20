#include <dice/diceMain.hxx> // Only include where main() function is, should be first

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "setup.hxx"

namespace glb = dice::glb; // shortened alias for global variable namespace

int main(int argc, char **argv)
{
  // Initialize the library
  dice::initialize(&argc, &argv);
  cfg::printDefines<dice::LOG_STD>(glb::console);

  // Add a few command line parameters
  int debugHook = glb::pParser->get<>("debugHook", dice::ParamsParser::defaultCategory(), 0,
                                      "Starts an infinite loop on startup to allow hooking a debugger");

  int debugWait = glb::pParser->get<>("debugWait", dice::ParamsParser::defaultCategory(), 0,
                                      "Wait for N seconds on startup to allow attaching a debugger");

  int help = glb::pParser->get<>("help", dice::ParamsParser::defaultCategory(), 0);

  int initOnly = glb::pParser->get<>("initOnly", dice::ParamsParser::defaultCategory(), 0,
                                     "Just initialize, report, and exit");

  int ignoreUnused = glb::pParser->get<>("ignoreUnusedParams", dice::ParamsParser::defaultCategory(), 0,
                                         "Do not exit when encountering unknown params");

  std::string restart = glb::pParser->get<>("restart", dice::ParamsParser::defaultCategory(), std::string(""),
                                            "The name of a restart file  (without the '.rst' or '_${RANK}.rst'(MPI) extension)");

  if ((initOnly) || (help))
    glb::console->setVerboseLevel(3);

  // Initialize the solver
  Problem *problem = new Problem(restart);

  // Report any parameter that was parsed but not used
  if (glb::pParser->reportUnused<dice::LOG_WARNING_SINGLE>())
  {
    // exit if required
    if (!ignoreUnused)
      return -1;
  }

  if (help)
    return 0;

  // Allow a debugger to hook
  if (debugHook)
    glb::mpiComWorld->debug();
  else if (debugWait > 0)
    glb::mpiComWorld->debug(debugWait);

  // Solve our problem !
  problem->solve();

  // Clean-up everything
  delete problem;

  dice::finalize();

  return 0;
}
