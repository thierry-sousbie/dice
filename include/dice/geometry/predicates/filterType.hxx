#ifndef __FILTER_TYPE_HXX__
#define __FILTER_TYPE_HXX__

/**
 * @file 
 * @brief  Filter types definition
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

namespace predicate { 
  /** \brief a namespace containing the definitions of the different predicate filter types */
  namespace filterType {

    /** \addtogroup Geometry
     *   \{
     */

    /** \enum Type Possible types of the predicate filters
     */
    enum Type {
      Raw=0,      /**< Fast inexact predicate, using only regular double precision */
      Exact=1,    /**< Slow exact predicate using GMP */
      Adaptive=2  /**< The best of both worlds ;) */
    };

    /** \}*/


  }
}

#include "../../internal/namespace.footer"
#endif
