#ifndef __INITIAL_CONDITIONS_INTERFACE_HXX__
#define __INITIAL_CONDITIONS_INTERFACE_HXX__

#include <dice/cosmo/cosmology.hxx>
#include <dice/mesh/tesselation/implicitTesselation.hxx>

#include "../units.hxx"

template <class M>
struct InitialConditionsInterfaceT
{
  typedef M Mesh;
  static const int NDIM=M::NDIM;
  static const int NDIM_W=M::NDIM_W;
  typedef dice::ImplicitTesselationT<NDIM,NDIM_W> ImplicitTesselation;
  
  typedef typename Mesh::Vertex Vertex;
  typedef typename Mesh::Simplex Simplex;
  
  struct Params
  {
    typedef dice::cosmo::Cosmology::Params CosmoParams;
    
    Params()
    {
      useCosmo=0;
      aStart=0;
      mass=1;

      G=6.67384e-11; // in m3 kg-1 s-2
      H=3.2409e-18;  // value in mks for h0=1 
      // NOTE: H is multiplied by h0 in the solver, except if H=1.

      x0.resize(NDIM_W,-1);
      delta.resize(NDIM_W,2);

      unitLength=1;
      unitVelocity=1;
      unitMass=1;    

      defaultUnitLength="1m";
      defaultUnitVelocity="1m/s";
      defaultUnitMass="1kg";
    }
  
    std::vector<double> x0;
    std::vector<double> delta;
    
    int useCosmo;
    double aStart;  
    double mass;
    double G;
    double H;
    CosmoParams cosmoParams;

    double unitLength;
    double unitVelocity;
    double unitMass;    

    std::string defaultUnitLength;
    std::string defaultUnitVelocity;
    std::string defaultUnitMass;
  };


  virtual std::string parserCategory() const {return "init";}
  
  virtual ~InitialConditionsInterfaceT()
  {}
  
  virtual ImplicitTesselation *createImplicitTesselation()=0;
  virtual void initialize(const double *x0,
			  const double *delta,
			  double t,
			  const Units &units, int pass)=0;  

  // Release temporary storage (this does not invalidate the reset the object!)
  virtual void release(){};
  //virtual void displace(double *coords)=0;
  virtual void displace(double *coords, double &density, Vertex *v)=0;
  virtual void displace(double *coords, Simplex *s)
  {
    double dummy;
    displace(coords,dummy,s->getVertex(0));
  }
 
  virtual long getNPasses() const {return 1;}
  virtual Params getParams() {return Params();}
  virtual bool useLagrangianMass() const {return false;}
  virtual std::string getName() const=0;
};

#endif
