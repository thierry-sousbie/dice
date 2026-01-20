#ifndef __BOOST_MULTIPRECISION_FOR_EIGEN_3_WRAPPER_HXX__
#define __BOOST_MULTIPRECISION_FOR_EIGEN_3_WRAPPER_HXX__

#include "cpp11SupportWrapper_nullptr.hxx" // eigen needs nullptr ...
#include <Eigen/Core>
#include <boost/multiprecision/number.hpp>

/**
 * @file
 * @brief This header should be included to enable usage of boost::multiprecision numbers
 * with eigen3 algebra library.
 * @author Thierry Sousbie
 */

namespace Eigen
{

  template <class Backend,
            boost::multiprecision::expression_template_option ExpressionTemplates>
  struct NumTraits<boost::multiprecision::number<Backend, ExpressionTemplates>>
      : GenericNumTraits<boost::multiprecision::number<Backend, ExpressionTemplates>>
  {
    typedef boost::multiprecision::number<Backend, ExpressionTemplates> MyNumberType;

    enum
    {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 10,
      AddCost = 10,
      MulCost = 40
    };

    typedef MyNumberType Real;
    typedef MyNumberType NonInteger;

    // Constants
    /*
    inline static Real Pi       (long Precision = MyNumberType::default_precision())
    {    return mpfr::const_pi(Precision);        }
    inline static Real Euler    (long Precision = MyNumberType::default_precision())
    {    return mpfr::const_euler(Precision);     }
    inline static Real Log2     (long Precision = MyNumberType::default_precision())
    {    return mpfr::const_log2(Precision);      }
    inline static Real Catalan  (long Precision = MyNumberType::default_precision())
    {    return mpfr::const_catalan(Precision);   }
    */
    // inline static Real epsilon  (long Precision = MyNumberType::default_precision())
    // {    return mpfr::machine_epsilon(Precision); }
    inline static Real epsilon(const Real &x)
    {
      return std::numeric_limits<MyNumberType>::epsilon();
    } // mpfr::machine_epsilon(x); }
    inline static Real epsilon()
    {
      return std::numeric_limits<MyNumberType>::epsilon();
    }

    inline static Real dummy_precision()
    {
      return epsilon() * 1.E4;
    }
  };

  namespace internal
  {
    /*
    template<> inline mpfr::mpreal random<mpfr::mpreal>()
    {
      return mpfr::random();
    }

    template<> inline mpfr::mpreal random<mpfr::mpreal>(const mpfr::mpreal& a, const mpfr::mpreal& b)
    {
      return a + (b-a) * random<mpfr::mpreal>();
    }

    inline bool isMuchSmallerThan(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& eps)
    {
      return mpfr::abs(a) <= mpfr::abs(b) * eps;
    }

    inline bool isApprox(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& eps)
    {
      return mpfr::isEqualFuzzy(a,b,eps);
    }

    inline bool isApproxOrLessThan(const mpfr::mpreal& a, const mpfr::mpreal& b, const mpfr::mpreal& eps)
    {
      return a <= b || mpfr::isEqualFuzzy(a,b,eps);
    }
    */

    /*
    template <class Backend,
        boost::multiprecision::expression_template_option ExpressionTemplates,
        class T>
    inline T cast< boost::multiprecision::number<Backend,ExpressionTemplates>,T>
    (const boost::multiprecision::number<Backend,ExpressionTemplates>& x)
    {
      return x.template convert_to<T>();
    }
    */
    /*
    template<> inline long double cast<mpfr::mpreal,long double>(const mpfr::mpreal& x)
    { return x.toLDouble(); }

    template<> inline double cast<mpfr::mpreal,double>(const mpfr::mpreal& x)
    { return x.toDouble(); }

    template<> inline long cast<mpfr::mpreal,long>(const mpfr::mpreal& x)
    { return x.toLong(); }

    template<> inline int cast<mpfr::mpreal,int>(const mpfr::mpreal& x)
    { return int(x.toLong()); }
    */

    /*

    // Specialize GEBP kernel and traits for mpreal (no need for peeling, nor complicated stuff)
    // This also permits to directly call mpfr's routines and avoid many temporaries produced by mpreal
    template<>
    class gebp_traits<mpfr::mpreal, mpfr::mpreal, false, false>
    {
    public:
      typedef mpfr::mpreal ResScalar;
      enum {
        nr = 1,
        mr = 1,
        LhsProgress = 1,
        RhsProgress = 1
      };
    };

    template<typename Index, bool ConjugateLhs, bool ConjugateRhs>
    struct gebp_kernel<mpfr::mpreal,mpfr::mpreal,Index,1,1,ConjugateLhs,ConjugateRhs>
    {
      typedef mpfr::mpreal mpreal;

      EIGEN_DONT_INLINE
      void operator()(mpreal* res, Index resStride, const mpreal* blockA, const mpreal* blockB, Index rows, Index depth, Index cols, mpreal alpha,
                      Index strideA=-1, Index strideB=-1, Index offsetA=0, Index offsetB=0)
      {
        if(rows==0 || cols==0 || depth==0)
          return;

        mpreal  acc1(0,mpfr_get_prec(blockA[0].mpfr_srcptr())),
    tmp (0,mpfr_get_prec(blockA[0].mpfr_srcptr()));

        if(strideA==-1) strideA = depth;
        if(strideB==-1) strideB = depth;

        for(Index i=0; i<rows; ++i)
    {
      for(Index j=0; j<cols; ++j)
        {
    mpreal *C1 = res + j*resStride;

    const mpreal *A = blockA + i*strideA + offsetA;
    const mpreal *B = blockB + j*strideB + offsetB;

    acc1 = 0;
    for(Index k=0; k<depth; k++)
      {
        mpfr_mul(tmp.mpfr_ptr(), A[k].mpfr_srcptr(), B[k].mpfr_srcptr(), mpreal::get_default_rnd());
        mpfr_add(acc1.mpfr_ptr(), acc1.mpfr_ptr(), tmp.mpfr_ptr(),  mpreal::get_default_rnd());
      }

    mpfr_mul(acc1.mpfr_ptr(), acc1.mpfr_srcptr(), alpha.mpfr_srcptr(), mpreal::get_default_rnd());
    mpfr_add(C1[i].mpfr_ptr(), C1[i].mpfr_srcptr(), acc1.mpfr_srcptr(),  mpreal::get_default_rnd());
        }
    }
      }
    };
    */
  } // end namespace internal
}

#endif
