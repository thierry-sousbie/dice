#ifndef __STATIC_POTENTIAL_SOLUTIONS_CHAOTIC_HXX__
#define __STATIC_POTENTIAL_SOLUTIONS_CHAOTIC_HXX__

#include "staticPotentialSolutionsInterface.hxx"

#include "../../internal/namespace.header" 

namespace cosmo {

  template <int ND>
  class StaticPotentialSolutions_ChaoticT : 
    public StaticPotentialSolutionsInterfaceT<ND>
  {
  public:
    typedef StaticPotentialSolutionsInterfaceT<ND> Base;

    static const int NDIM=ND;

    template <class R, class PM>
    StaticPotentialSolutions_ChaoticT(R *reader, PM &paramsManager, double mass, double G_)
    {
      initialize(reader,paramsManager,mass,G_);
    }
  
    double acc_over_r(int i, double r, double *coords) const
    {
      double x=coords[0];
      double y=coords[1];
      double x2=x*x;
      double y2=y*y;
      //double r2=r*r;
      double result;

      if (i==0)
	{
	  double up=-q2*(-2*Re*r+3*x2+y2);
	  double down=2*r*(Re*Rc2*q2-q2*x2*r+q2*y2*r+Re*q2*x2+Re*y2);
	  result=up/down;
	}
      else
	{
	  double up=(q2*x2+3*q2*y2+2*Re*r);
	  double down=2*r*(Re*Rc2*q2-q2*x2*r+q2*y2*r+Re*q2*x2+Re*y2);
	  result=up/down;
	}

      if (r==0)
	{
	  if (i==0)
	    result = -1/Rc2;
	  else
	    result = 1.0/(Rc2*q2);
	}

      return -result;
      /*
      double result=1.0/((2*r*Rc2*q2*Re+Re*x2*q2+Re*y2-pow(r,3)*q2)*r);
      
      if (i==0)
	result*=x*(2*r*Re-3*x2-y2)*q2;
      else
	result*=y*(2*r*Re+q2*x2+3*q2*y2);

      if (r==0)
	{
	  if (i==0)
	    result = pMass*(2*Re*q2)/(2*Rc2*q2*Re);
	  else
	    result = pMass*(2*Re)/(2*Rc2*q2*Re);
	}

      return -result;
      */
    }

    double density(double r) const
    {
      // double x2=coords[0]*coords[0];
      // double y2=coords[1]*coords[1];
      // return 0.5*log(Rc2+x2+y2*q2_inv-r*(x2-y2)*epsilon);
      return 0;
    }

    template <class R, class PM>
    void initialize(R *reader, PM &paramsManager,double mass, double G_)
    {
      Rc = 0.2;
      Rc = paramsManager.
	get("Rc",Base::parserCategory(),Rc,reader,
	    PM::FILE_FIRST,
	    "Central radius");

      q = 0.9;
      q = paramsManager.
	get("q",Base::parserCategory(),q,reader,
	    PM::FILE_FIRST,
	    "Excentricity ?");

      epsilon = 0.5;
      epsilon = paramsManager.
	get("epsilon",Base::parserCategory(),epsilon,reader,
	    PM::FILE_FIRST,
	    "Perturbation amplitude (epsilon=0 -> regular orbits)");
   
      initialize(mass,G_);
    }

  private:

    void initialize(double mass, double G_)
    {
      //const double pi = 4.0*atan(1.0);

      G=G_;
      pMass=mass;    

      Rc2=Rc*Rc;
      q2=q*q;
      q2_inv=1.0/q2;
      if (epsilon==0) Re=1.E100;
      else Re=1.0/epsilon;      
    }

    double pMass;
    double Rc,q,epsilon;    
    double Rc2,q2,q2_inv;    
    double Re;
    double G;   
  };
}

#include "../../internal/namespace.footer"
#endif
