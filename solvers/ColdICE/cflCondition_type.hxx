#ifndef __CFL_CONDITION_TYPE_HXX__
#define __CFL_CONDITION_TYPE_HXX__

#include <dice/tools/types/typeSelect.hxx>

struct CflConditionTypeV {
  enum Type{NONE=0, 
	    RHOMAX=(1<<0), 
	    CFL=(1<<1), 
	    CFL_RHOMAX=(1<<0)|(1<<1), 
	    UNDEFINED=100};
};

typedef CflConditionTypeV::Type CflConditionType;

struct CflConditionTypeSelect : public dice::TypeSelectT<CflConditionTypeV> {
  CflConditionTypeSelect()
  {
    this->insert("NONE",CflConditionTypeV::NONE);
    this->insert("RHOMAX",CflConditionTypeV::RHOMAX);    
    this->insert("CFL_RHOMAX",CflConditionTypeV::CFL_RHOMAX);    
    this->insert("CFL",CflConditionTypeV::CFL);    
  }
  std::string name() {return "cfl_condition_type";}
};

#endif
