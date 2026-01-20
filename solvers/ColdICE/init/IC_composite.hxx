#ifndef __IC_COMPOSITE_HXX__
#define __IC_COMPOSITE_HXX__

#include <math.h>

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include <dice/mesh/tesselation/compositeTesselation.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_composite : public InitialConditionsInterfaceT<M>
{
private:
  static const int NDIM = M::NDIM;
  static const int NDIM_W = M::NDIM_W;

public:
  typedef InitialConditionsInterfaceT<M> Interface;
  typedef typename Interface::ImplicitTesselation ImplicitTesselation;
  typedef typename Interface::Params Params;
  typedef typename ImplicitTesselation::Cell Cell;
  typedef typename Interface::Vertex Vertex;

  static std::string name() { return "composite"; }
  std::string getName() const { return name(); }

  template <class PM, class R>
  IC_composite(PM &paramsManager, R *reader,
               dice::MpiCommunication *com = 0) :                         // nullptr
                                                  compositeTesselation(0) // nullptr
  {
    for (int i = 0; i < NDIM; ++i)
    {
      x0[i] = -1;
      delta[i] = 2;
    }
    for (int i = NDIM; i < NDIM_W; ++i)
    {
      x0[i] = 0;
      delta[i] = 0;
    }
    parse(paramsManager, reader);
  }

  void initialize(const double *x0, const double *delta, double t,
                  const Units &units, int pass)
  {
    if (pass == 0)
    {
      release();
      compositeTesselation = new CompositeTesselation(x0, delta);
      for (unsigned long i = 0; i < tesselations.size(); ++i)
        compositeTesselation->add(tesselations[i]);
    }
  }

  void release()
  {
    if (compositeTesselation != 0) // nullptr
    {
      delete compositeTesselation;
      compositeTesselation = 0; // nullptr
    }
  }

  void displace(double *coords, double &density, Vertex *v)
  {
    typename CompositeTesselation::Cell cell(0, v->getGeneration().id());
    long component = compositeTesselation->getComponentId(cell);

    density = relativeDensity[component];

    Components &vel = velVector[component];
    for (int i = NDIM; i < NDIM_W; ++i)
      coords[i] += vel.val[i - NDIM];

    if (angMomentum[component] != 0)
    {
      Components &c = center[component];

      double v[NDIM];
      if (NDIM == 2)
      {
        v[0] = -(coords[1] - c.val[1]) * angMomentum[component];
        v[1] = (coords[0] - c.val[0]) * angMomentum[component];
      }
      else
      {
        // Components &mom=angMomentumVec[component];
      }

      for (int i = 0; i < NDIM; ++i)
        coords[i + NDIM] += v[i];
    }
  }

  Params getParams()
  {
    Params p;

    std::copy_n(x0, NDIM_W, &p.x0.front());
    std::copy_n(delta, NDIM_W, &p.delta.front());

    p.useCosmo = 0;
    p.cosmoParams.set(1.0, 0, 0, 1.0);
    p.mass = 1;
    p.aStart = 1.e-2;

    p.unitLength = 1.0;
    p.unitVelocity = 1.0;
    p.unitMass = 1.0;

    p.defaultUnitLength = "-";
    p.defaultUnitVelocity = "-";
    p.defaultUnitMass = "-";

    p.G = 1;
    p.H = 1;

    return p;
  }

  ImplicitTesselation *createImplicitTesselation()
  {
    CompositeTesselation *result = new CompositeTesselation(x0, delta);

    for (unsigned long i = 0; i < tesselations.size(); ++i)
      result->add(tesselations[i]);

    return static_cast<ImplicitTesselation *>(result);
  }

private:
  typedef dice::SimplicialGridT<NDIM, NDIM_W, Cell, false> SimplicialGrid;
  typedef typename SimplicialGrid::Params SimplicialGridParams;
  typedef dice::CompositeTesselationT<NDIM, NDIM_W, Cell> CompositeTesselation;

  CompositeTesselation *compositeTesselation;
  std::vector<ImplicitTesselation *> tesselations;
  double x0[NDIM_W];
  double delta[NDIM_W];

  struct Components
  {
    double val[NDIM];
    void normalize()
    {
      double tmp = val[0] * val[0];
      for (int i = 1; i < NDIM; ++i)
        tmp += val[i] * val[i];
      if (tmp != 0)
      {
        tmp = sqrt(tmp);
        for (int i = 0; i < NDIM; ++i)
          val[i] /= tmp;
      }
    }
  };

  std::vector<Components> velVector;
  std::vector<Components> angMomentumVec;
  std::vector<double> angMomentum;
  std::vector<Components> center;
  std::vector<double> relativeDensity;

  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {
    std::string newComponent = "none";

    // Using "bBox.." to avoid parameter conflicts
    for (int i = 0; i < NDIM_W; i++)
    {
      x0[i] = paramsManager.get("bBoxX0", Interface::parserCategory(), x0[i], i,
                                reader, PM::FILE_FIRST,
                                "Initial coordinates of the bounding box lower left corner");

      delta[i] = paramsManager.get("bBoxDelta", Interface::parserCategory(), delta[i], i,
                                   reader, PM::FILE_FIRST,
                                   "Size of the bounding box");
    }

    int pass = 0;
    int nBox = 0;
    do
    {
      newComponent = "none";
      newComponent = paramsManager.get("add", Interface::parserCategory(), newComponent, pass,
                                       reader, PM::FILE_FIRST,
                                       "The type of component to add ['box'].");

      if (newComponent != "none")
      {
        velVector.resize(velVector.size() + 1);
        angMomentumVec.resize(angMomentum.size() + 1);

        velVector[pass].val[0] = (pass > 0) ? velVector[pass - 1].val[0] : 0;
        velVector[pass].val[0] = paramsManager.get("velVector", Interface::parserCategory(),
                                                   velVector[pass].val[0], pass * NDIM + 0,
                                                   reader, PM::FILE_FIRST,
                                                   "Global velocity vector of the component.");

        relativeDensity.push_back((pass > 0) ? relativeDensity[pass - 1] : 1.0);
        relativeDensity[pass] = paramsManager.get("relativeDensity", Interface::parserCategory(),
                                                  relativeDensity[pass], pass,
                                                  reader, PM::FILE_FIRST,
                                                  "Relative density of the component (w.r.t other components). The value is only meaningfull w.r.t. other components value as the total mass is normalized (see solver.mass).");

        for (int i = 1; i < NDIM; ++i)
        {
          velVector[pass].val[i] =
              (pass > 0) ? velVector[pass - 1].val[i] : velVector[pass].val[0];
          velVector[pass].val[i] = paramsManager.get("velVector", Interface::parserCategory(),
                                                     velVector[pass].val[i], pass * NDIM + i,
                                                     reader, PM::FILE_FIRST,
                                                     "Global velocity vector of the component.");
        }

        if (NDIM > 2)
        {
          angMomentumVec[pass].val[0] = (pass > 0) ? angMomentumVec[pass - 1].val[0] : 0;
          angMomentumVec[pass].val[0] = paramsManager.get("angularMomentumVec", Interface::parserCategory(),
                                                          angMomentumVec[pass].val[0], pass * NDIM + 0,
                                                          reader, PM::FILE_FIRST,
                                                          "Angular momentum vector of the component.");

          for (int i = 1; i < NDIM; ++i)
          {
            angMomentumVec[pass].val[i] =
                (pass > 0) ? angMomentumVec[pass - 1].val[i] : angMomentumVec[pass].val[0];
            angMomentumVec[pass].val[i] = paramsManager.get("angularMomentumVec", Interface::parserCategory(),
                                                            angMomentumVec[pass].val[i], pass * NDIM + i,
                                                            reader, PM::FILE_FIRST,
                                                            "Angular momentum vector of the component.");
          }
          angMomentumVec[pass].normalize();
        }

        angMomentum.push_back((pass > 0) ? angMomentum[pass - 1] : 0);
        angMomentum[pass] = paramsManager.get("angularMomentum", Interface::parserCategory(),
                                              angMomentum[pass], pass,
                                              reader, PM::FILE_FIRST,
                                              "Angular momentum amplitude of the component.");
      }

      center.resize(center.size() + 1);
      if (newComponent == "box")
      {
        SimplicialGridParams sgParams;
        sgParams.parse(reader, paramsManager, Interface::parserCategory(), nBox++);
        tesselations.push_back(new SimplicialGrid(sgParams));

        for (int i = 0; i < NDIM; ++i)
          center[pass].val[i] = sgParams.x0[i] + sgParams.delta[i] * 0.5;
      }
      else if (newComponent != "none")
      {
        dice::PRINT_SRC_INFO(dice::LOG_ERROR);
        dice::glb::console->print<dice::LOG_ERROR>("Invalid component type '%s'.\n", newComponent.c_str());
        throw std::runtime_error("Invalid parameter specified when creating 'composite' initial conditions.\n");
      }
      ++pass;
    } while (newComponent != "none");
  }
};

#endif
