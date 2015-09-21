#ifndef LIBQD_QD_BACKEND_FOR_BOOST_MULTI_HXX__
#define LIBQD_QD_BACKEND_FOR_BOOST_MULTI_HXX__

#include <limits>
#include <math.h>

#ifdef HAVE_QD
#include <boost/multiprecision/detail/float_string_cvt.hpp>
#include <qd/qd_real.h>

namespace boost{
  namespace multiprecision{
    namespace backends{
      
      struct qd_backend
      {
	typedef mpl::list<signed char, short, int, long, long long>   signed_types;
	typedef mpl::list<unsigned char, unsigned short, 
			  unsigned int, unsigned long, unsigned long long> unsigned_types;
	typedef mpl::list<float, double, long double> float_types;
	typedef int exponent_type;

      private:
	qd_real m_value;
      public:
	qd_backend() BOOST_NOEXCEPT : m_value(0) {}
	qd_backend(const qd_backend& o) BOOST_NOEXCEPT : m_value(o.m_value) {}
	qd_backend& operator = (const qd_backend& o) BOOST_NOEXCEPT
	{
	  m_value = o.m_value;
	  return *this;
	}
	template <class T>
	qd_backend(const T& i, const typename enable_if_c<is_convertible<T, qd_real>::value>::type* = 0) BOOST_NOEXCEPT
	  : m_value(i) {}
	template <class T>
	typename enable_if_c<is_arithmetic<T>::value || is_convertible<T, qd_real>::value, qd_backend&>::type operator = (const T& i) BOOST_NOEXCEPT
	{
	  m_value = i;
	  return *this;
	}
	qd_backend& operator = (const char* s)
	{
	  boost::multiprecision::detail::convert_from_string(*this, s);
	  return *this;
	}
	void swap(qd_backend& o) BOOST_NOEXCEPT
	{
	  std::swap(m_value, o.value());
	}
	std::string str(std::streamsize digits, std::ios_base::fmtflags f)const
	{
	  return boost::multiprecision::detail::convert_to_string
	    (*this, digits ? digits : 37, f);
	}
	void negate() BOOST_NOEXCEPT
	{
	  m_value = -m_value;
	}
	int compare(const qd_backend& o)const
	{
	  return m_value == o.m_value ? 0 : m_value < o.m_value ? -1 : 1;
	}
	template <class T>
	int compare(const T& i)const
	{
	  return m_value == i ? 0 : m_value < i ? -1 : 1;
	}
	qd_real& value()
	{
	  return m_value;
	}
	const qd_real& value()const
	{
	  return m_value;
	}
      };

      inline void eval_add(qd_backend& result, const qd_backend& a)
      {
	result.value() += a.value();
      }
      template <class A>
      inline void eval_add(qd_backend& result, const A& a)
      {
	result.value() += a;
      }
      inline void eval_subtract(qd_backend& result, const qd_backend& a)
      {
	result.value() -= a.value();
      }
      template <class A>
      inline void eval_subtract(qd_backend& result, const A& a)
      {
	result.value() -= a;
      }
      inline void eval_multiply(qd_backend& result, const qd_backend& a)
      {
	result.value() *= a.value();
      }
      template <class A>
      inline void eval_multiply(qd_backend& result, const A& a)
      {
	result.value() *= a;
      }
      inline void eval_divide(qd_backend& result, const qd_backend& a)
      {
	result.value() /= a.value();
      }
      template <class A>
      inline void eval_divide(qd_backend& result, const A& a)
      {
	result.value() /= a;
      }

      inline void eval_add(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = a.value() + b.value();
      }
      template <class A>
      inline void eval_add(qd_backend& result, const qd_backend& a, const A& b)
      {
	result.value() = a.value() + b;
      }
      inline void eval_subtract(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = a.value() - b.value();
      }
      template <class A>
      inline void eval_subtract(qd_backend& result, const qd_backend& a, const A& b)
      {
	result.value() = a.value() - b;
      }
      template <class A>
      inline void eval_subtract(qd_backend& result, const A& a, const qd_backend& b)
      {
	result.value() = a - b.value();
      }
      inline void eval_multiply(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = a.value() * b.value();
      }
      template <class A>
      inline void eval_multiply(qd_backend& result, const qd_backend& a, const A& b)
      {
	result.value() = a.value() * b;
      }
      inline void eval_divide(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = a.value() / b.value();
      }

      template <class R>
      inline void eval_convert_to(R* result, const qd_backend& val)
      {
	*result = 
	  static_cast<R>(val.value().x[0])+static_cast<R>(val.value().x[1])+
	  static_cast<R>(val.value().x[2])+static_cast<R>(val.value().x[3]);
      }

      inline void eval_frexp(qd_backend& result, const qd_backend& arg, int* exp)
      {	
	::frexp(arg.value().x[0], exp);
	result.value()=mul_pwr2(arg.value(),::ldexp(1.0,-(*exp)));
      }

      inline void eval_ldexp(qd_backend& result, const qd_backend& arg, int exp)
      {
	result.value() = mul_pwr2(arg.value(),::ldexp(1.0,exp));
      }

      inline void eval_floor(qd_backend& result, const qd_backend& arg)
      {
	result.value() = floor(arg.value());
      }
      inline void eval_ceil(qd_backend& result, const qd_backend& arg)
      {
	result.value() = ceil(arg.value());
      }
      inline void eval_sqrt(qd_backend& result, const qd_backend& arg)
      {
	result.value() = sqrt(arg.value());
      }
      inline int eval_fpclassify(const qd_backend& arg)
      {
	return isnan(arg.value()) ? FP_NAN : isinf(arg.value()) ? FP_INFINITE : arg.value() == 0 ? FP_ZERO : FP_NORMAL;
      }

      inline void eval_increment(qd_backend& arg)
      {
	arg.value()+=1;
      }
      inline void eval_decrement(qd_backend& arg)
      {
	arg.value()-=1;
      }

      /*********************************************************************
       *
       * abs/fabs:
       *
       *********************************************************************/

      inline void eval_abs(qd_backend& result, const qd_backend& arg)
      {
	result.value() = fabs(arg.value());
      }
      inline void eval_fabs(qd_backend& result, const qd_backend& arg)
      {
	result.value() = fabs(arg.value());
      }

      /*********************************************************************
       *
       * Floating point functions:
       *
       *********************************************************************/
      /*
	inline void eval_trunc(qd_backend& result, const qd_backend& arg)
	{
	if(isnanq(arg.value()) || isinfq(arg.value()))
	{
	result = boost::math::policies::raise_rounding_error
	("boost::multiprecision::trunc<%1%>(%1%)", 0, 
	number<qd_backend, et_off>(arg), 
	number<qd_backend, et_off>(arg), 
	boost::math::policies::policy<>()).backend();
	return;
	}
	result.value() = trunc(arg.value());
	}
      */

      /*
      // 
      // This doesn't actually work... rely on our own default version instead.
      //
      inline void eval_round(qd_backend& result, const qd_backend& arg)
      {
      if(isnanq(arg.value()) || isinf(arg.value()))
      {
      result = boost::math::policies::raise_rounding_error(
      "boost::multiprecision::trunc<%1%>(%1%)", 0, 
      number<qd_backend, et_off>(arg), 
      number<qd_backend, et_off>(arg), 
      boost::math::policies::policy<>()).backend();
      return;
      }
      result.value() = roundq(arg.value());
      }
      */

      inline void eval_exp(qd_backend& result, const qd_backend& arg)
      {
	result.value() = exp(arg.value());
      }
      inline void eval_log(qd_backend& result, const qd_backend& arg)
      {
	result.value() = log(arg.value());
      }
      inline void eval_log10(qd_backend& result, const qd_backend& arg)
      {
	result.value() = log10(arg.value());
      }
      inline void eval_sin(qd_backend& result, const qd_backend& arg)
      {
	result.value() = sin(arg.value());
      }
      inline void eval_cos(qd_backend& result, const qd_backend& arg)
      {
	result.value() = cos(arg.value());
      }
      inline void eval_tan(qd_backend& result, const qd_backend& arg)
      {
	result.value() = tan(arg.value());
      }
      inline void eval_asin(qd_backend& result, const qd_backend& arg)
      {
	result.value() = asin(arg.value());
      }
      inline void eval_acos(qd_backend& result, const qd_backend& arg)
      {
	result.value() = acos(arg.value());
      }
      inline void eval_atan(qd_backend& result, const qd_backend& arg)
      {
	result.value() = atan(arg.value());
      }
      inline void eval_sinh(qd_backend& result, const qd_backend& arg)
      {
	result.value() = sinh(arg.value());
      }
      inline void eval_cosh(qd_backend& result, const qd_backend& arg)
      {
	result.value() = cosh(arg.value());
      }
      inline void eval_tanh(qd_backend& result, const qd_backend& arg)
      {
	result.value() = tanh(arg.value());
      }
      /*
	inline void eval_fmod(qd_backend& result, const qd_backend& a, const qd_backend& b)
	{
	result.value() = fmod(a.value(), b.value());
	}
      */
      inline void eval_pow(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = pow(a.value(), b.value());
      }
      inline void eval_atan2(qd_backend& result, const qd_backend& a, const qd_backend& b)
      {
	result.value() = atan2(a.value(), b.value());
      }

      inline bool eval_is_zero(const qd_backend& val) BOOST_NOEXCEPT
      {      
	return boost::multiprecision::default_ops::eval_is_zero(val);
      }

      inline int eval_get_sign(const qd_backend& val) BOOST_NOEXCEPT
      {      
	return boost::multiprecision::default_ops::eval_get_sign(val);
      }

    } // namespace backends    
  } // namespaces
} // namespaces

namespace std{
  template <boost::multiprecision::expression_template_option ExpressionTemplates>
  class numeric_limits<boost::multiprecision::number<boost::multiprecision::backends::qd_backend, ExpressionTemplates> > : public numeric_limits<qd_real>
  {};
}

#endif

#endif

