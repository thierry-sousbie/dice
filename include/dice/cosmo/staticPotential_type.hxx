#ifndef __STATIC_POTENTIAL_TYPE_HXX__
#define __STATIC_POTENTIAL_TYPE_HXX__

#include "../tools/types/typeSelect.hxx"

#include "../internal/namespace.header" 

namespace cosmo {

  struct StaticPotentialTypeV {
    enum Type{PLUMMER=0,UDF=1,CHAOTIC=2,UNDEFINED=100};
  };

  typedef StaticPotentialTypeV::Type StaticPotentialType;

  struct StaticPotentialTypeSelect : 
    public TypeSelectT<StaticPotentialTypeV> {
    StaticPotentialTypeSelect()
    {
      this->insert("PLUMMER",StaticPotentialTypeV::PLUMMER);
      this->insert("UDF",StaticPotentialTypeV::UDF);   
      this->insert("CHAOTIC",StaticPotentialTypeV::CHAOTIC);   
    }
    std::string name() {return "staticPotentialType";}
  };

} // namespace cosmo

#include "../internal/namespace.footer"
#endif
