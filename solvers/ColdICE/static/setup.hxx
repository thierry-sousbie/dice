#ifndef __COLDICE_SETUP_HXX__
#define __COLDICE_SETUP_HXX__

#include <dice/solver/solverInterface.hxx>

#include "config.h"
#include "../coldice.hxx"

#include "staticPotential.hxx"

// boundary conditions (can be PERIODIC or BOXED)
#ifdef D_BOUNDARY_TYPE
#undef D_BOUNDARY_TYPE
#define D_BOUNDARY_TYPE BOXED
#endif

typedef 
Coldice< D_DIMS_COUNT , dice::BoundaryType:: D_BOUNDARY_TYPE>
SolverImpl;


typedef dice::slv::SolverInterfaceT< StaticPotential<SolverImpl> > Problem;

#endif
