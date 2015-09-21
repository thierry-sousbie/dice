#ifndef __IC_SINE_WAVES_HXX__
#define __IC_SINE_WAVES_HXX__

#include <math.h>

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_sineWaves : 
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

  static std::string name() {return "sineWaves";}
  std::string getName() const {return name();}

  template <class PM, class R>
  IC_sineWaves(PM &paramsManager, R *reader, 
	       dice::MpiCommunication *com=0) //nullptr
  {
    parse(paramsManager,reader);
  }

  void initialize(const double *x0, const double *delta, double t, 
		  const Units &units, int pass)
  {
    std::copy_n(x0,NDIM,x0_);
    std::copy_n(delta,NDIM,delta_);

    vFact[0]=0;
    vFact[1]=0;

    xFact[0]=1.0;
    xFact[1]=0.0;

    if (cosmo)
      {
	double aStart=units.cosmology.a_of_tau(t,0);	

	vFact[0]=
	  units.cosmology.h_over_h0(aStart)*
	  units.cosmology.f1(aStart)*
	  units.H*units.length*aStart/units.velocity; 	

	vFact[1]=
	  units.cosmology.h_over_h0(aStart)*
	  units.cosmology.f2(aStart)*
	  units.H*units.length*aStart/units.velocity; 

	xFact[0]=units.cosmology.d1_plus(aStart);
	xFact[1]=units.cosmology.d2_plus(aStart);
      }    
  }

  void displace(double *coords, double &density, Vertex *v)
  {
    double x[NDIM];
    for (int i=0;i<NDIM;++i) 
      x[i]=(coords[i]-x0_[i])/delta_[i];
    
    for (int i=0;i<NDIM;++i) 
      {	     
	double dx = xFact[0]*dx1(x,i) + xFact[1]*dx2(x,i);
	double dv = xFact[0]*vFact[0]*dx1(x,i) + xFact[1]*vFact[1]*dx2(x,i);

	coords[i]+=dx*delta_[i];
	coords[i+NDIM]+=dv*delta_[i];
      }
     
  }

  Params getParams()
  {
    Params p;
    std::copy_n(sgParams.x0,NDIM,&p.x0.front());
    std::copy_n(sgParams.delta,NDIM,&p.delta.front());

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
    
    return p;
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
    cosmo = 1;
    cosmo = paramsManager.
      get("cosmo",Interface::parserCategory(),cosmo,
	  reader,PM::FILE_FIRST,
	  "Whether to use cosmology");

    if (cosmo)
      {
	std::fill_n(sgParams.x0,NDIM,-0.5);
	std::fill_n(sgParams.delta,NDIM,1.0);
      }
    else
      {
	std::fill_n(sgParams.x0,NDIM,-0.5);
	std::fill_n(sgParams.delta,NDIM,1.0);
      }
    sgParams.parse(reader,paramsManager,Interface::parserCategory());

    amplitude[0] = (cosmo)?2:0.1;
    amplitude[0] = paramsManager.
      get("amplitude",Interface::parserCategory(),amplitude[0],0,
	  reader,PM::FILE_FIRST,
	  "Amplitude of the sine waves");

    for (int i=1;i<NDIM;++i)
      {
	amplitude[i] = amplitude[i-1];
	amplitude[i] = paramsManager.
	  get("amplitude",Interface::parserCategory(),amplitude[i],i,
	      reader,PM::FILE_FIRST,
	      "Amplitude of the sine waves");
      }   
  }

  double dx1(double x[NDIM], int index)
  {
    static const double twoPi = 2.0*M_PI;    
    return -sin(2.0*M_PI*(x[index]-0.5))*amplitude[index]/twoPi;
  }

  double dx2(double x[NDIM], int index)
  {
    // Removed second order !
    return 0;

    static const double twoPi = 2.0*M_PI;
    
    if (!cosmo) return 0;

    double result = 1.0/(4.0*twoPi)*sin(2.0L*M_PI*(x[index]-0.5))*amplitude[index];

    double tmp=0;
    for (int i=0;i<NDIM;++i)
      if (i!=index) tmp+=cos(2.0L*M_PI*(x[i]-0.5))*amplitude[i];
    
    return result*tmp;
  }

  int cosmo;
  double amplitude[NDIM];
  double x0_[NDIM];
  double delta_[NDIM]; 
  double vFact[2];
  double xFact[2];
};

#endif
