#ifndef __SIMPLEX_TYPE_HXX__
#define __SIMPLEX_TYPE_HXX__

/**
 * @file 
 * @brief Defines Simplex type implementation specialization template parameters
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

/**
 * \namespace simplexType
 * \brief Defines Simplex type implementation specialization template parameters
 */
namespace simplexType {

enum Type {
  VerticesOnly=0,
  WithSegments=1
};

};

/** \}*/
#include "../internal/namespace.footer"
#endif
