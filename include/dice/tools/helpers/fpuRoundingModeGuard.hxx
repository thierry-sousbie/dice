#ifndef FPU_ROUNDING_MODE_GUARD_HXX__
#define FPU_ROUNDING_MODE_GUARD_HXX__

/**
 * @file
 * @brief Defines a class that needs to be overloaded to set appropriate fpu rounding
 * depending on the data type.
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */

namespace hlp
{
  template <class T>
  class FpuRoundingModeGuardT
  {
  public:
    // Avoid warnings
    FpuRoundingModeGuardT()
    {
    }
    ~FpuRoundingModeGuardT()
    {
    }
  };
}

/** \}*/
#include "../../internal/namespace.footer"
#endif
