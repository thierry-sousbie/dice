#ifndef INITIAL_CONDITIONS_HXX__
#define INITIAL_CONDITIONS_HXX__

#include <dice/tools/helpers/typeList.hxx>

#include "initialConditionsFactory.hxx"
#include "initialConditionsInterface.hxx"

// Include initial conditions headers
#include "IC_uniformGrid.hxx"
#include "IC_sineWaves.hxx"
#include "IC_phasedWave.hxx"
#include "IC_cosmo.hxx"
#include "IC_composite.hxx"
// #include "IC_dummy.hxx"

// Register initial conditions, add your own conditions to this list !
template <class M>
class InitialConditionsListT : public dice::hlp::TypeListT<

                                   IC_uniformGrid<M>,
                                   IC_sineWaves<M>,
                                   IC_phasedWave<M>,
                                   IC_cosmo<M>,
                                   IC_composite<M>
                                   //,IC_dummy<M>

                                   >
{
};

#endif
