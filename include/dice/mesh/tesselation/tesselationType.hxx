#ifndef __TESSELATION_TYPE_HXX__
#define __TESSELATION_TYPE_HXX__

#include "../../tools/types/typeSelect.hxx"
/**
 * @file
 * @brief Defines an enum and its map to strings for the different types of tesselation of
 * a regular grid into simplices
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

struct TesselationTypeV
{
  enum Type
  {
    ANY = 0,
    MINIMAL = 1,
    REGULAR = 2,
    ALTERNATE = 3,
    USER_DEFINED = 4,
    UNDEFINED = 100
  };
};

typedef TesselationTypeV::Type TesselationType;

struct TesselationTypeSelect : public TypeSelectT<TesselationTypeV>
{
  TesselationTypeSelect()
  {
    this->insert("REGULAR", TesselationTypeV::REGULAR);
    this->insert("ALTERNATE", TesselationTypeV::ALTERNATE);
    this->insert("MINIMAL", TesselationTypeV::MINIMAL);
    this->insert("USER_DEFINED", TesselationTypeV::USER_DEFINED);
    this->insert("ANY", TesselationTypeV::ANY);
  }
  std::string name() { return "tesselationType"; }
};

/** \}*/
#include "../../internal/namespace.footer"
#endif
