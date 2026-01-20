#ifndef IC_DUMMY_HXX__
#define IC_DUMMY_HXX__

// This file is a dummy initial condition implementation designed to be used as a template.
// Copy/paste it to a IC_<your_type>.hxx file, replace all instances of "dummy"
// with <your_type> and fill the gaps as needed. See other IC_[...].hxx files to check
// how they are implementated. You can then edit initialConditions.hxx to add your
// own IC type.

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_dummy : public InitialConditionsInterfaceT<M>
{
private:
  static const int NDIM = M::NDIM;
  static const int NDIM_W = M::NDIM_W;

public:
  typedef InitialConditionsInterfaceT<M> Interface;
  typedef typename Interface::ImplicitTesselation ImplicitTesselation;
  typedef typename Interface::Params Params;
  typedef typename Interface::Vertex Vertex;

  static std::string name() { return "dummy"; }
  std::string getName() const { return name(); }

  template <class PM, class R>
  IC_dummy(PM &paramsManager, R *reader,
           dice::MpiCommunication *com = 0)
  {
    parse(paramsManager, reader);
  }

  // This will be called before displace is ever called.
  // Use it to initialize any data you'll need to displace the vertices
  void initialize(const double *x0, const double *delta, double t,
                  const Units &units, int pass)
  {
    // Copy bounding box information to use in displace
    std::copy_n(x0, NDIM, this->x0);
    std::copy_n(delta, NDIM, this->delta);
  }

  // Implement the initial conditions here
  // -> On each call, update the new coordinates (coords) and
  // density (density) associated to the current vertex (v)
  // Input coordinates are the ones given by the implicit tesselation
  // you provided in createImplicitTesselation(), and input density is
  // set to the average density.
  // Note that density will be multiplied by the projected or lagrangian
  // volume of the simplices afterward to compute the masses depending on
  // useLagrangianMass() function (see below).
  // See e.g. IC_sineWaves.hxx, IC_cosmo.hxx, ...
  void displace(double *coords, double &density, Vertex *v)
  {
    // Transform coords to be in the range [0,1]
    double x[NDIM];
    for (int i = 0; i < NDIM; ++i)
      x[i] = (coords[i] - x0_[i]) / delta_[i];

    // set the new coordinates here, and update density if needed
    /*
    for (int i=0;i<NDIM;++i)
      {
  coords[i] = ... ;
      }
    */
  }

  // Here the initial conditions default paramters can be set if needed
  // Default values will be provided automatically if you do not do it here ...
  // See the implementation of Params in initialConditionsInterface.hxx
  Params getParams()
  {
    Params p;
    // Copy the the user specified bounding box of the implicit tesselation
    // to IC parameters
    std::copy_n(sgParams.x0, NDIM, &p.x0.front());
    std::copy_n(sgParams.delta, NDIM, &p.delta.front());

    // Exemple of default parameters equivalent to the ones which
    // are automatically provided
    /*
    p.useCosmo=cosmo;
    p.cosmoParams.set(1.0,0,0,1.0);
    p.mass=1;
    p.aStart=1.e-2;

    p.unitLength=1.0;
    p.unitVelocity=1.0;
    p.unitMass=1.0;

    p.defaultUnitLength="-";
    p.defaultUnitVelocity="-";
    p.defaultUnitMass="-";

    p.G=1;
    p.H=1;
    */
    return p;
  }

  // This should return an implicit tesselation as is.
  // You probably do not need to modifiy this
  ImplicitTesselation *createImplicitTesselation()
  {
    SimplicialGrid *sg = new SimplicialGrid(sgParams);
    return reinterpret_cast<ImplicitTesselation *>(sg);
  }

  // If true, initial mass of simplex S is computed as
  //  M_S=(M_tot/V_tot * lagrangian_volume(S))
  // If false, then
  //  M_S=(M_tot/V_tot * projected_volume(S))
  bool useLagrangianMass() const { return true; }

  // Set this to the number of passes that should be used to
  // displace vertices
  // See e.g. IC_cosmo.hxx
  long getNPasses() const { return 1; }

  // Here you can release memory as needed once displacing is complete ...
  // Note that the destructor will be called correctly anyway, but memory is
  // just released earlier if you do it here.
  // Don't bother unless you are using a lot of memory !
  // See e.g. IC_cosmo.hxx
  void release()
  {
  }

private:
  typedef typename ImplicitTesselation::Cell Cell;

  // This is the default type of grid we will generate.
  // In this case, periodicity is the same as the mesh, but you can set the last
  // template parameter to dice::BoundaryType::PERIODIC or dice::BoundaryType::BOXED
  typedef dice::SimplicialGridT<NDIM, NDIM_W, Cell, M::IS_PERIODIC> SimplicialGrid;
  typedef typename SimplicialGrid::Params SimplicialGridParams;

  // Parameters of the implicit grid
  SimplicialGridParams sgParams;

  // Here we do the parsing of required parameters
  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {
    // This will parse default parameters for the implicit grid
    sgParams.parse(reader, paramsManager, Interface::parserCategory());

    // Add other parser parameters here as needed
    // This parameter for instance can be used to query whether cosmological
    // conditions shoudl be used
    /*
    cosmo = 1;
    cosmo = paramsManager.
      get("cosmo",Interface::parserCategory(),cosmo,
    reader,PM::FILE_FIRST,
    "Whether to use cosmology");
    */
  }

  // int cosmo;
  double x0[NDIM];
  double delta[NDIM];
};

#endif
