#ifndef __REGULAR_GRID_SYMMETRY_ENUM_HXX__
#define __REGULAR_GRID_SYMMETRY_ENUM_HXX__

#include "../tools/types/typeSelect.hxx"

/**
 * @file 
 * @brief Defines an enum for different grid symmetries
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

struct RGSV {
  enum Type{
    NONE      = (0),    

    PLANAR0   = (1<<0),
    PLANAR1   = (1<<1),
    PLANAR01  = (1<<1)|(1<<0),
    PLANAR2   = (1<<2),
    PLANAR02  = (1<<2)|(1<<0),
    PLANAR12  = (1<<2)|(1<<1),   
    PLANAR012 = (1<<2)|(1<<1)|(1<<0),

    AXIAL0    = (1<<3),
    AXIAL1    = (1<<4),
    AXIAL01   = (1<<4)|(1<<3),
    AXIAL2    = (1<<5),
    AXIAL02   = (1<<5)|(1<<3),
    AXIAL12   = (1<<5)|(1<<4),   
    AXIAL012  = (1<<5)|(1<<4)|(1<<3), 
   
    CONSTANT0 = (1<<6),
    CONSTANT1 = (1<<7),
    CONSTANT2 = (1<<8),

    CENTRAL   = (1<<9),

    UNDEFINED= NONE
  };  
};

typedef typename RGSV::Type RegularGridSymmetryE;

struct RegularGridSymmetrySelect : public TypeSelectT<RGSV> {
  RegularGridSymmetrySelect()
  {
    this->insert("NONE",RGSV::NONE);

    this->insert("PLANAR0",RGSV::PLANAR0);
    this->insert("PLANAR1",RGSV::PLANAR1);
    this->insert("PLANAR2",RGSV::PLANAR2);
    this->insert("PLANAR01",RGSV::PLANAR01);
    this->insert("PLANAR02",RGSV::PLANAR02);
    this->insert("PLANAR12",RGSV::PLANAR12);
    this->insert("PLANAR012",RGSV::PLANAR012);
    
    this->insert("AXIAL0",RGSV::AXIAL0);
    this->insert("AXIAL1",RGSV::AXIAL1);
    this->insert("AXIAL2",RGSV::AXIAL2);
    this->insert("AXIAL01",RGSV::AXIAL01);
    this->insert("AXIAL02",RGSV::AXIAL02);
    this->insert("AXIAL12",RGSV::AXIAL12);
    this->insert("AXIAL012",RGSV::AXIAL012);

    this->insert("CONSTANT0",RGSV::CONSTANT0);
    this->insert("CONSTANT1",RGSV::CONSTANT1);
    this->insert("CONSTANT2",RGSV::CONSTANT2);

    this->insert("CENTRAL",RGSV::CENTRAL);
  }
  std::string name() {return "RegularGridSymmetry";}
};

/** \}*/
#include "../internal/namespace.footer"
#endif
