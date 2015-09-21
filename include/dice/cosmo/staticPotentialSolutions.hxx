#ifndef __STATIC_POTENTIAL_SOLUTIONS_HXX__
#define __STATIC_POTENTIAL_SOLUTIONS_HXX__

#include "staticPotential_type.hxx"
#include "staticPotentials/staticPotentialSolutions_UDF.hxx"
#include "staticPotentials/staticPotentialSolutions_Plummer.hxx"
#include "staticPotentials/staticPotentialSolutions_Chaotic.hxx"

#include "../internal/namespace.header" 

namespace cosmo {

  template <int ND>
  class StaticPotentialSolutionsT
  {
    
  public:
    typedef StaticPotentialType  PotentialType;
    typedef StaticPotentialTypeV PotentialTypeV;
    typedef StaticPotentialTypeSelect PotentialTypeSelect;

    static const int NDIM = ND;
    /*
    template <class R, class PM>
    StaticPotentialSolutionsT(StaticPotentialType t=StaticPotentialTypeV::PLUMMER, 
			      double mass=1, double G_=1)
    {
      interface=NULL;
      switch (t)
	{
	case StaticPotentialTypeV::UDF:
	  interface = new StaticPotentialSolutions_UDFT<ND>(reader,manager,mass,G_);
	case StaticPotentialTypeV::PLUMMER:
	  interface = new StaticPotentialSolutions_PlummerT<ND>(reader,manager,mass,G_);
	default:
	  {};
	}      
    }
    */

    StaticPotentialSolutionsT()
    {}

    template <class R, class PM>
    StaticPotentialSolutionsT(StaticPotentialType t,
			      R *reader, PM &paramsManager,
			      double mass=1, double G_=1)
    {
      selectPotential(t,reader,paramsManager,mass,G_);
    }
    
    template <class R, class PM>
    void selectPotential(StaticPotentialType t,
			 R *reader, PM &manager,
			 double mass=1, double G_=1)
    {
      if (interface!=NULL) delete interface;
      interface=NULL;
      switch (t)
	{
	case StaticPotentialTypeV::UDF:
	  interface = new StaticPotentialSolutions_UDFT<ND>(reader,manager,mass,G_);
	  break;
	case StaticPotentialTypeV::PLUMMER:
	  interface = new StaticPotentialSolutions_PlummerT<ND>(reader,manager,mass,G_);
	  break;
	case StaticPotentialTypeV::CHAOTIC:
	  interface = new StaticPotentialSolutions_ChaoticT<ND>(reader,manager,mass,G_);
	  break;
	default:
	  {};
	}      
    }

    ~StaticPotentialSolutionsT()
    {
      delete interface;
    }

    
    double acc_over_r(int i, double r, double *coords) const
    {
      return interface->acc_over_r(i,r,coords);
    }

    double density(double r) const
    {
      return interface->density(r);
    }
 
  private:

    StaticPotentialSolutionsInterfaceT<ND> *interface;
  };
}

#include "../internal/namespace.footer"
#endif
