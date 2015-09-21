#ifndef __MEMORY_POOL_ITERATOR_POLICIES_HXX__
#define __MEMORY_POOL_ITERATOR_POLICIES_HXX__

/**
 * @file 
 * @brief Defines polcies for the memory pool iterators
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 *   \{
 */

/** 
 * \namespace iteratorThreadModel 
 * \brief memoryPool iterators threading policies
 */
namespace iteratorThreadModel {
  class Blocks{}; /*!< each thread is given a contiguous block of elements */
  class Alternate{}; /*!< each thread is given one element every num_threads elements */
  class SortedBlocks{}; /*!< each thread is given two blocks (each of them contiguous): one in the sorted region, one in the unsorted region  */
}

/** \}*/
#include "../../internal/namespace.footer"
#endif
