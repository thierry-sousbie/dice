#ifndef __BOUNDARY_TYPE_HXX__
#define __BOUNDARY_TYPE_HXX__

/**
 * @file
 * @brief An enum for boundary conditions types
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup Geometry
 *   \{
 */

/** \namespace BoundaryType
 * \brief a namespace containing the definitions of the different boundary types
 */
namespace BoundaryType
{

  /** \enum Type possible types of boundary conditions */
  enum Type
  {
    NONE = 0,    /**< No boundary conditions (i.e. unbounded) */
    BOXED = 1,   /**< space is limited by a bounding box */
    PERIODIC = 2 /**< Periodic boundary conditions */
  };

};

/** \}*/
#include "../internal/namespace.footer"
#endif
