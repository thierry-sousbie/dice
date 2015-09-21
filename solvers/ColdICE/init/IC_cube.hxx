#ifndef __IC_CUBE_HXX__
#define __IC_CUBE_HXX__

#include <math.h>

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_cube : 
  public InitialConditionsInterfaceT<M>
{
private:
  static const int NDIM=M::NDIM;
  static const int NDIM_W=M::NDIM_W;

public:
  typedef InitialConditionsInterfaceT<M> Interface;
  typedef typename Interface::ImplicitTesselation ImplicitTesselation;
  typedef typename Interface::Params Params;
  typedef typename Interface::Vertex Vertex;
  
  static std::string name() {return "cube";}
  std::string getName() const {return name();}

  template <class PM, class R>
  IC_cube(PM &paramsManager, R *reader, 
	  dice::MpiCommunication *com=0) //nullptr
  {
    parse(paramsManager,reader);
  }

  void initialize(const double *x0,const double *delta, double t, 
		  const Units &units, int pass)
  {
    std::copy_n(x0,NDIM,x0_);
    std::copy_n(delta,NDIM,delta_);    
  }

  void displace(double *coords, double &density, Vertex *v)
  {
    double x[NDIM];
    for (int i=0;i<NDIM;++i) 
      x[i]=(coords[i]-x0_[i])/delta_[i];

    for (int i=0;i<NDIM;++i) 
      {
	coords[i]=(x[i]*scale_[i]) + center_[i];
	coords[i+NDIM]=velocity_[i];
      }
  }

  ImplicitTesselation *createImplicitTesselation()
  { 
    SimplicialGrid *sg = new SimplicialGrid(sgParams);
    return reinterpret_cast<ImplicitTesselation*>(sg);
  }

  bool useLagrangianMass() const {return true;}

private:
  typedef typename ImplicitTesselation::Cell Cell;
  typedef dice::SimplicialGridT<NDIM,NDIM_W,Cell,M::IS_PERIODIC> SimplicialGrid;
  typedef typename SimplicialGrid::Params SimplicialGridParams;
  
  SimplicialGridParams sgParams;

  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {      
    sgParams.parse(reader,paramsManager,Interface::parserCategory());
    
    scale_[0] = 0.2;
    scale_[0] = paramsManager.
      get("scale",Interface::parserCategory(),scale_[0],0,
	  reader,PM::FILE_FIRST,
	  "Scale of the cube (w.r.t. the domain size)");

    for (int i=1;i<NDIM;++i)
      {
	scale_[i] = scale_[0];
	scale_[i] = paramsManager.
	  get("scale",Interface::parserCategory(),scale_[i],i,
	      reader,PM::FILE_FIRST,
	      "Scale of the cube (w.r.t. the domain size)");
      }

    center_[0] = 0.2;
    center_[0] = paramsManager.
      get("center",Interface::parserCategory(),center_[0],0,
	  reader,PM::FILE_FIRST,
	  "Center of the sine waves");

    for (int i=1;i<NDIM;++i)
      {
	center_[i] = center_[0];
	center_[i] = paramsManager.
	  get("center",Interface::parserCategory(),center_[i],i,
	      reader,PM::FILE_FIRST,
	      "Position of the cube");
      }

    velocity_[0] = 0.2;
    velocity_[0] = paramsManager.
      get("velocity",Interface::parserCategory(),velocity_[0],0,
	  reader,PM::FILE_FIRST,
	  "Velocity vector");

    for (int i=1;i<NDIM;++i)
      {
	velocity_[i] = velocity_[0];
	velocity_[i] = paramsManager.
	  get("velocity",Interface::parserCategory(),velocity_[i],i,
	      reader,PM::FILE_FIRST,
	      "Velocity of the cube");
      }
  }
  
  double scale_[3];
  double center_[3];
  double velocity_[3];
  double x0_[3];
  double delta_[3];  
};

#endif
