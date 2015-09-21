#ifndef BOOST_MULTIPRECISION_QD_WRAPPER_HXX__
#define BOOST_MULTIPRECISION_QD_WRAPPER_HXX__

/**
 * @file 
 * @brief This header defines a DoubleDouble and QuadDouble data types using 
 * boost::multiprecision and double and quadruple double rpecision numbers from libqd.
 * @author Thierry Sousbie
 */

#include "internal/qdBackend.hxx"
#include "internal/ddBackend.hxx"
#include "../helpers/fpuRoundingModeGuard.hxx"

namespace boost{
  namespace multiprecision{
    namespace backends{

      struct dd_backend;
      struct qd_backend;

    }
    using backends::dd_backend;
    using backends::qd_backend;

    template<>
    struct number_category<backends::dd_backend> : 
      public mpl::int_<number_kind_floating_point> {};

    template<>
    struct number_category<dd_real> : 
      public mpl::int_<number_kind_floating_point> {};

    typedef number<dd_backend, et_off> doubledouble;

    template<>
    struct number_category<backends::qd_backend> : 
      public mpl::int_<number_kind_floating_point> {};

    template<>
    struct number_category<qd_real> : 
      public mpl::int_<number_kind_floating_point> {};

    typedef number<qd_backend, et_off> quaddouble;
  }
}

#include "../../internal/namespace.header"
typedef boost::multiprecision::doubledouble DoubleDouble;
typedef boost::multiprecision::quaddouble   QuadDouble;

namespace hlp {
  template <>
  class FpuRoundingModeGuardT<DoubleDouble>
  {
  public:
    FpuRoundingModeGuardT()
    {fpu_fix_start(&old_cw);}
  
    ~FpuRoundingModeGuardT()
    {fpu_fix_end(&old_cw);}
  private:
    unsigned int old_cw;
  };


  template <>
  class FpuRoundingModeGuardT<QuadDouble>
  {
  public:
    FpuRoundingModeGuardT()
    {fpu_fix_start(&old_cw);}
    
    ~FpuRoundingModeGuardT()
    {fpu_fix_end(&old_cw);}
  private:
    unsigned int old_cw;
  };
}

#include "../../internal/namespace.footer"
#endif
