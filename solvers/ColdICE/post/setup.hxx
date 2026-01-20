#ifndef __COLDICE_SETUP_HXX__
#define __COLDICE_SETUP_HXX__

#include <dice/solver/solverInterface.hxx>

#include "../config.h"
#include "../coldice.hxx"

typedef Coldice<D_DIMS_COUNT, dice::BoundaryType::D_BOUNDARY_TYPE>
    SolverImpl;

#include "post.hxx"
typedef dice::slv::SolverInterfaceT<ColdicePost<SolverImpl>> Problem;

#endif
