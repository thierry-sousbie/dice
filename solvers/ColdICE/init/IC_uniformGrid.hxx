#ifndef __IC_UNIFORM_GRID_HXX__
#define __IC_UNIFORM_GRID_HXX__

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_uniformGrid : public InitialConditionsInterfaceT<M>
{
private:
  static const int NDIM = M::NDIM;
  static const int NDIM_W = M::NDIM_W;

public:
  typedef InitialConditionsInterfaceT<M> Interface;
  typedef typename Interface::ImplicitTesselation ImplicitTesselation;
  typedef typename Interface::Params Params;
  typedef typename Interface::Vertex Vertex;

  static std::string name() { return "uniformGrid"; }
  std::string getName() const { return name(); }

  template <class PM, class R>
  IC_uniformGrid(PM &paramsManager, R *reader,
                 dice::MpiCommunication *com = 0) // nullptr
  {
    parse(paramsManager, reader);
  }

  void initialize(const double *x0, const double *delta, double t,
                  const Units &units, int pass)
  {
  }

  void displace(double *coords, double &density, Vertex *v)
  {
  }

  Params getParams()
  {
    Params p;
    std::copy_n(sgParams.x0, NDIM, &p.x0.front());
    std::copy_n(sgParams.delta, NDIM, &p.delta.front());
    return p;
  }

  ImplicitTesselation *createImplicitTesselation()
  {
    SimplicialGrid *sg = new SimplicialGrid(sgParams);
    return reinterpret_cast<ImplicitTesselation *>(sg);
  }

  bool useLagrangianMass() const { return true; }

private:
  typedef typename ImplicitTesselation::Cell Cell;
  typedef dice::SimplicialGridT<NDIM, NDIM_W, Cell, M::IS_PERIODIC> SimplicialGrid;
  typedef typename SimplicialGrid::Params SimplicialGridParams;

  SimplicialGridParams sgParams;

  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {
    sgParams.parse(reader, paramsManager, Interface::parserCategory());
  }
};

#endif
