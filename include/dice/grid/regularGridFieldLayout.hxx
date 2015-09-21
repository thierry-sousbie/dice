#ifndef __REGULAR_GRID_FIELD_LAYOUT_HXX__
#define __REGULAR_GRID_FIELD_LAYOUT_HXX__

/**
 * @file 
 * @brief An enum use to define fields memory layout for regular grids.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
   *   \{
   */

/** \namespace RegularGridFieldLayout
 *   \brief a namespace containing definitions for regular grid fields memory layouts 
*/
namespace regularGridFieldLayout {
  
  /** \enum Type possible types of regular grid fields memory layouts */
  enum Type {
    CONSECUTIVE = 0, /**< each field is stored as a consecutive array (X1X2..Y1Y2..Z1Z2..)*/
    INTERLEAVED = 1  /**< the nth elements of each field are grouped (X1Y1Z1X2Y2Z2...) */
  };

};
/** \}*/

#include "../internal/namespace.footer"
#endif
