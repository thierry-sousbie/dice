#ifndef __GRID_KERNELS_OPERATIONS_HXX__
#define __GRID_KERNELS_OPERATIONS_HXX__

/**
 * @file
 * @brief Definition possible operations over grid kernels
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

namespace gridKernel
{

  enum KernelOp
  {
    IMPRINT = 0, /**< Set pixels overlaped by the kernel to a given value */
    APPLY = 1,   /**< Apply the kernel */
    COPY = 2,    /**< Copy the value overlapped by the kernel */
    INTERNAL1    /**< internal operation, do not use ;) */
  };

}

/** \}*/
#include "../internal/namespace.footer"
#endif
