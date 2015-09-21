#ifndef __STATIC_POTENTIAL_SOLUTIONS_INTERFACE_HXX__
#define __STATIC_POTENTIAL_SOLUTIONS_INTERFACE_HXX__

#include "../../internal/namespace.header" 

namespace cosmo {

  template <int ND>
  class StaticPotentialSolutionsInterfaceT
  {
  public:
    static std::string parserCategory() {return "potential";}

    virtual ~StaticPotentialSolutionsInterfaceT()
    {}

    virtual double acc_over_r(int i, double r, double *coords) const=0;
    //virtual double density(double r) const=0;
  };
  
}

#include "../../internal/namespace.footer"
#endif
