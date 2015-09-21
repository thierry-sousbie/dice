#ifndef __VAL_LOCATION_TYPE_HXX__
#define __VAL_LOCATION_TYPE_HXX__

#include "../tools/types/typeSelect.hxx"

#include "../internal/namespace.header"

struct ValLocationTypeV {
  enum Type {CELL=0, VERTEX=1, UNDEFINED=-1};
};

typedef ValLocationTypeV::Type ValLocationType;

struct ValLocationTypeSelect : public TypeSelectT<ValLocationTypeV> {
  ValLocationTypeSelect()
  {
    insert("cell",ValLocationTypeV::CELL);
    insert("vertex",ValLocationTypeV::VERTEX);
  }
  std::string name() {return "val_location_type";}
};

#include "../internal/namespace.footer"
#endif
