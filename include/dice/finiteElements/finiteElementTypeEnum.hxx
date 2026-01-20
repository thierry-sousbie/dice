#ifndef __FINITE_ELEMENT_TYPE_ENUM_HXX__
#define __FINITE_ELEMENT_TYPE_ENUM_HXX__

/**
 * @file
 * @brief  Define an enum for finite elements  types
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup FiniteElements
 *   \{
 */

//! protects the definition of finite element type as an enum
namespace fETypeE
{
  /// defines different types of finite element
  enum Type
  {
    simplex,
    cuboid
  };
}

/** \}*/
#include "../internal/namespace.footer"
#endif
