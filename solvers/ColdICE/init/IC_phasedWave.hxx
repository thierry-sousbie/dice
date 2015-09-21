#ifndef __IC_PHASED_WAVE_HXX__
#define __IC_PHASED_WAVE_HXX__

#include <math.h>

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_phasedWave : 
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

  static std::string name() {return "phasedWave";}
  std::string getName() const {return name();}

  template <class PM, class R>
  IC_phasedWave(PM &paramsManager, R *reader, 
		dice::MpiCommunication *com=0) //nullptr
  {
    parse(paramsManager,reader);
  }

  void initialize(const double *x0, const double *delta, double t,
		  const Units &units, int pass)
  {
    std::copy_n(x0,NDIM,x0_);
    std::copy_n(delta,NDIM,delta_);    
    
    aStart = units.cosmology.a_of_tau(t,0);
    
    double waveAmp = 1.0/aCross;
    for (int i=0;i<NDIM;++i)
      initAmp[i] = waveAmp*units.cosmology.d1_plus(aStart)*delta[i];

    //vFact = 100.0*pow(aStart,-1.5)*sqrt(aStart);    // for oM=1
    vFact=
      units.cosmology.h_over_h0(aStart)*  
      units.cosmology.f1(aStart)*    
      units.H*units.length*aStart/units.velocity; 
    
    //printf("vFact = %lg / %lg\n",vFact,100.0*pow(aStart,-1.5)*sqrt(aStart));
  }

  void displace(double *coords, double &density, Vertex *v)
  {    
    static const double twopi = 2.0*M_PI;
    
    double x[NDIM];
    for (int i=0;i<NDIM;++i) 
      x[i]=(coords[i]-x0_[i])/delta_[i];

    double tmp=twopi*kp*x[0];
    for (int i=0;i<NDIM-1;++i)
       tmp+=epsilon[i]*(kp*kp)/(k[i]*k[i])*cos(twopi*k[i]*x[i+1]);
    tmp=sin(tmp);

    coords[0+NDIM]=initAmp[0]/(twopi*kp)*tmp;
    for (int i=1;i<NDIM;++i)
      coords[i+NDIM]=-initAmp[i]/(twopi*k[i-1])*epsilon[i-1]*sin(twopi*k[i-1]*x[i])*tmp;
    
    for (int i=0;i<NDIM;++i)
      {
	coords[i]+=coords[i+NDIM];
	coords[i+NDIM]*=vFact;
	//printf("V->%lg\n",coords[i+NDIM]);	
      }
    
    /*
    double vz = initAmp/(twopi*kp)*
      sin(twopi*kp*qz + epsilon_a*(kp*kp)/(ka*ka)*cos(twopi*ka*qy));
    double vy = -initAmp/(twopi*ka)*epsilon_a*
      sin(twopi*kp*qz + epsilon_a*(kp*kp)/(ka*ka)*cos(twopi*ka*qy))*
      sin(twopi*ka*qy);
    double vx = 0.0;

    velocity[3*ii+0] = vx*vfact;
    velocity[3*ii+1] = vy*vfact;
    velocity[3*ii+2] = vz*vfact;
				
    position[3*ii+0] = boxlength*qx + vx;
    position[3*ii+1] = boxlength*qy + vy;
    position[3*ii+2] = boxlength*qz + vz;
    */
  }

  Params getParams()
  {
    Params p;

    std::copy_n(sgParams.x0,NDIM,&p.x0.front());
    std::copy_n(sgParams.delta,NDIM,&p.delta.front());
    
    p.useCosmo=1;
    p.cosmoParams.set(1.0,0,0,1.0);
    p.mass=1;
    p.aStart=1.0/(1.0+19.0);

    p.unitLength=3.085678e22;    // 1.0 Mpc
    p.unitVelocity=1.e3;         // 1km/s
    p.unitMass=1.989e40;         // 1e10 solar masses

    p.defaultUnitLength="1Mpc";
    p.defaultUnitVelocity="1km/s";
    p.defaultUnitMass="1E10 M_sol";
    
    p.G=6.67384e-11; // in m3 kg-1 s-2
    p.H=3.2409e-18;  // value in mks for h0=1

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

    std::fill_n(sgParams.x0,NDIM,0.0);
    std::fill_n(sgParams.delta,NDIM,100.0);
    
    sgParams.parse(reader,paramsManager,Interface::parserCategory());
    
    aCross = 0.15;
    aCross = paramsManager.
      get("aCross",Interface::parserCategory(),aCross,
	  reader,PM::FILE_FIRST,
	  "Time of first shell crossing");

    kp = 1.0;
    kp = paramsManager.
      get("kp",Interface::parserCategory(),kp,
	  reader,PM::FILE_FIRST,
	  "Plane wave parameter");

    k[0] = 1.0;
    k[0] = paramsManager.
      get("k",Interface::parserCategory(),k[0],0,
	  reader,PM::FILE_FIRST,
	  "Perturbation parameter");

    epsilon[0] = 0.5;
    epsilon[0] = paramsManager.
      get("epsilon",Interface::parserCategory(),epsilon[0],0,
	  reader,PM::FILE_FIRST,
	  "");
    
    for (int i=1;i<NDIM-1;++i)
      {
	k[i] = k[i-1];
	k[i] = paramsManager.
	  get("k",Interface::parserCategory(),k[i],i,
	      reader,PM::FILE_FIRST,
	      "Perturbation parameter");

	epsilon[i] = epsilon[i-1]+0.2; 
	epsilon[i] = paramsManager.
	  get("epsilon",Interface::parserCategory(),epsilon[i],i,
	      reader,PM::FILE_FIRST,
	      "");
      }
  }
 
  double k[NDIM-1];       // perturbation
  double epsilon[NDIM-1]; // 
  double kp;         // plane wave
  double vFact;

  double aStart;  
  double aCross;
  double initAmp[NDIM];
  double x0_[3];
  double delta_[3];  
};

#endif
