#ifndef __HELPERS_HXX__
#define __HELPERS_HXX__

#include <algorithm>
#include <cmath>
#include "../wrappers/cpp11SupportWrapper_utility.hxx"
// #include <boost/multiprecision/float128.hpp>

/**
 * @file
 * @brief Defines various helper functions that help a lot, especially when CPP11 is not
 * available ...
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 * \{
 */

/**
 * \namespace hlp
 * \brief Contains definitions of various helper functions/classes that help a lot !
 */

namespace hlp
{

  // **** SWAP_FLAG
  template <typename TF>
  void swapFlag(TF &val, TF flag)
  {
    if (val & flag)
      val &= (~flag);
    else
      val |= flag;
  }

  template <typename TF>
  void swapFlagValues(TF &val, TF flag1, TF flag2)
  {
    if ((val & flag1) == (val & flag2))
      return;
    hlp::swapFlag(val, flag1);
    hlp::swapFlag(val, flag2);
  }

  // **** FIRST_FIRST_DIFFERENCE
  template <class T, int N>
  struct FindFirstDifference
  {
    static int find(T A[N], T B[N])
    {
      int j;
      for (int i = 0; i < N; i++)
      {
        for (j = 0; j < N; j++)
        {
          if (A[i] == B[j])
            break;
        }
        if (j == N)
          return i;
      }
      return -1;
    }
  };

  template <class T>
  struct FindFirstDifference<T, 4>
  {
    static int find(T A[4], T B[4])
    {
      if ((A[0] != B[0]) && (A[0] != B[1]) && (A[0] != B[2]) && (A[0] != B[3]))
        return 0;
      if ((A[1] != B[0]) && (A[1] != B[1]) && (A[1] != B[2]) && (A[1] != B[3]))
        return 1;
      if ((A[2] != B[0]) && (A[2] != B[1]) && (A[2] != B[2]) && (A[2] != B[3]))
        return 2;
      if ((A[3] != B[0]) && (A[3] != B[1]) && (A[3] != B[2]) && (A[3] != B[3]))
        return 3;
      return -1;
    }
  };
  template <class T>
  struct FindFirstDifference<T, 3>
  {
    static int find(T A[3], T B[3])
    {
      if ((A[0] != B[0]) && (A[0] != B[1]) && (A[0] != B[2]))
        return 0;
      if ((A[1] != B[0]) && (A[1] != B[1]) && (A[1] != B[2]))
        return 1;
      if ((A[2] != B[0]) && (A[2] != B[1]) && (A[2] != B[2]))
        return 2;
      return -1;
    }
  };
  template <class T>
  struct FindFirstDifference<T, 2>
  {
    static int find(T A[2], T B[2])
    {
      if ((A[0] != B[0]) && (A[0] != B[1]))
        return 0;
      if ((A[1] != B[0]) && (A[1] != B[1]))
        return 1;
      return -1;
    }
  };
  template <class T>
  struct FindFirstDifference<T, 1>
  {
    static int find(T A[1], T B[1])
    {
      if (A[0] != B[0])
        return 0;
      return -1;
    }
  };

  // **** COMPARE_POINTERS_TO

  template <class C>
  struct comparePointersToLess : public std::binary_function<C *, C *, bool>
  {
    bool operator()(const C *A, const C *B) const
    {
      return (*A) < (*B);
    }
  };

  template <class C>
  struct comparePointersToMore : public std::binary_function<C *, C *, bool>
  {
    bool operator()(const C *A, const C *B) const
    {
      return (*A) > (*B);
    }
  };

  // **** CONSTANT_VALUE

  template <long N>
  struct ConstantValue
  {
    // static const long value = N;
    enum
    {
      value = N
    };
  };

  template <typename T>
  struct ConstantType
  {
    typedef T Type;
  };

  // Factorial

  template <long N>
  struct FactorialT
  {
    enum
    {
      value = N * FactorialT<N - 1>::value
    };
  };

  template <>
  struct FactorialT<0>
  {
    enum
    {
      value = 1
    };
  };

  // True / False types

  struct IsEnabled
  {
    static bool isEnabled() { return true; }
    enum
    {
      value = 1
    };
  };

  struct IsDisabled
  {
    static bool isEnabled() { return false; }
    enum
    {
      value = 0
    };
  };

  struct IsTrue
  {
    static bool isTrue() { return true; }
    enum
    {
      value = 1
    };
  };

  struct IsFalse
  {
    static bool isTrue() { return false; }
    enum
    {
      value = 0
    };
  };

  template <bool b>
  struct IsTrueT
  {
    typedef IsTrue Result;
  };

  template <>
  struct IsTrueT<false>
  {
    typedef IsFalse Result;
  };

  // make a class or type const without the ugliness of const_cast or having to define a new variable
  template <class T>
  static const T &makeConst(T &w) { return w; }

  // Type comparaison (N.B. : a type and its const version are considered the same)
  template <typename T, typename U>
  struct SameType
  {
    typedef typename SameType<const T, const U>::Result Result;
    enum
    {
      value = Result::value
    };
  };

  template <typename T, typename U>
  struct SameType<const T, const U>
  {
    typedef IsFalse Result;
    enum
    {
      value = Result::value
    };
  };

  template <typename T>
  struct SameType<const T, const T>
  {
    typedef IsTrue Result;
    enum
    {
      value = Result::value
    };
  };

  // IF_
  template <bool condition, class A, class B>
  struct IF_
  {
    typedef A Result;
  };

  template <class A, class B>
  struct IF_<false, A, B>
  {
    typedef B Result;
  };

  // MIN / MAX
  template <long A, long B>
  struct MinT
  {
    enum
    {
      value =
          IF_<(A <= B), ConstantValue<A>, ConstantValue<B>>::Result::value
    };
  };

  template <long A, long B>
  struct MaxT
  {
    enum
    {
      value =
          IF_<(A > B), ConstantValue<A>, ConstantValue<B>>::Result::value
    };
  };

  // INTEGER power
  template <long A, long B, long value>
  struct IntPowerHelper
  {
    typedef typename IntPowerHelper<A, B - 1, value * A>::Result Result;
  };

  template <long A, long value>
  struct IntPowerHelper<A, 0, value>
  {
    typedef ConstantValue<value> Result;
  };

  template <long A, long B>
  struct IntPower
  {
    typedef typename IntPowerHelper<A, B, 1>::Result Result;
    enum
    {
      value = Result::value
    };
  };

  // Get the sign of any type
  template <typename T>
  T sgn(T val)
  {
    return (T(0) < val) - (val < T(0));
  }

  // T is a std::vector
  template <typename T>
  struct IsVector
  {
    typedef IsFalse Result;
  };
  template <typename T, typename A>
  struct IsVector<std::vector<T, A>>
  {
    typedef IsTrue Result;
  };
  template <typename T, typename A>
  struct IsVector<const std::vector<T, A>>
  {
    typedef IsTrue Result;
  };

  // check if T is a primitive type (i.e. char, int, float ...)
  template <class T>
  struct IsPrimitiveType
  {
    typedef IsFalse Result;
  };

  template <class T>
  struct IsPrimitiveType<T *>
  {
    typedef IsTrue Result;
  };

  template <class T>
  struct IsPrimitiveType<T &>
  {
    typedef typename IsPrimitiveType<T>::Result Result;
  };

  template <class T>
  struct IsPrimitiveType<const T>
  {
    typedef typename IsPrimitiveType<T>::Result Result;
  };

  template <class T>
  struct IsPrimitiveType<const T *>
  {
    typedef IsTrue Result;
  };

  template <class T>
  struct IsPrimitiveType<const T &>
  {
    typedef typename IsPrimitiveType<T>::Result Result;
  };

  template <>
  struct IsPrimitiveType<bool>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<char>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<short>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<int>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<long>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<long long>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<unsigned char>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<unsigned short>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<unsigned int>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<unsigned long>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<unsigned long long>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<float>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<double>
  {
    typedef IsTrue Result;
  };
  template <>
  struct IsPrimitiveType<long double>
  {
    typedef IsTrue Result;
  };

  // BoostMultiStaticCast: cast a boost::multiprecision or a regular type with the same
  // interface
  namespace internal
  {
    // compact version of SFINAE (test if result_type exists)
    template <class T>
    struct IsBoostMultiExprTemplate
    {
      // boost::multi expression templates define a U::result_type
      template <class U>
      static char (&test(typename U::result_type const *))[1];
      // capture anything but U::result_type.
      template <class U>
      static char (&test(...))[2];

      // value is 1 if T::result_type is defined
      static const bool value = (sizeof(test<T>(0)) == 1);
    };

    template <class T1, class T2, bool is_float>
    struct MostPreciseFP_helper;

    template <class T1, class T2>
    struct MostPreciseFP_helper<T1, T2, false>
    {
      // Compile time error message
      struct ERROR_TYPE_MUST_BE_FLOAT;
      typedef ERROR_TYPE_MUST_BE_FLOAT Result;
    };

    template <class T1, class T2>
    struct MostPreciseFP_helper<T1, T2, true>
    {
      typedef typename IF_<
          (std::numeric_limits<T1>::digits > std::numeric_limits<T2>::digits),
          T1, T2>::Result Result;
    };

    // Just a regular static cast
    template <class OT, class I, bool P, bool ET>
    struct NumericStaticCastT
    {
      static OT cast(const I &val)
      {
        return static_cast<OT>(val);
      }
    };

    // boost::multi (multiprecision number)
    template <class OT, class I>
    struct NumericStaticCastT<OT, I, false, false>
    {

      static OT cast(const I &val)
      {
        return static_cast<OT>(val);
      }

      /*
  // This is needed when compiling with C++03
      static OT cast(const I &val)
      {return val.template convert_to<OT>();}
      */
    };

    // boost::multi (expression template)
    template <class OT, class I>
    struct NumericStaticCastT<OT, I, false, true>
    {
      static OT cast(const I &val)
      {
        typedef typename I::result_type Result;
        Result tmp = val; // evaluate the expression
        return NumericStaticCastT<OT, Result, false, false>::cast(tmp);
      }
    };

  }

  // Used to test if a type T is defined (n.b. it must be declared)
  template <typename T>
  struct HasDestructor
  {
    /* Has destructor :) */
    template <typename A>
    static IsTrue test(decltype(cpp11::declval<A>().~A()) *)
    {
      return IsTrue();
    }

    /* Has no destructor :( */
    template <typename A>
    static IsFalse test(...)
    {
      return IsFalse();
    }

    /* This will be either `std::true_type` or `std::false_type` */
    typedef decltype(test<T>(0)) type;

    static const bool value = type::value; /* Which is it? */
  };

  // Get the floating point type with highest precision
  template <class T1, class T2>
  struct MostPreciseFP
  {
    typedef typename internal::MostPreciseFP_helper<
        T1,
        T2,
        (!std::numeric_limits<T1>::is_integer) &&
            (!std::numeric_limits<T2>::is_integer)>::Result Result;
  };

  template <class OT, class I>
  OT numericStaticCast(const I &val)
  {
    typedef typename IsPrimitiveType<I>::Result IsPrimitive;
    typedef typename internal::IsBoostMultiExprTemplate<I> IsExpr;
    typedef internal::NumericStaticCastT<OT, I, IsPrimitive::value, IsExpr::value>
        NumericStaticCast;

    return NumericStaticCast::cast(val);
  }

  // Check if a T is empty, using empty base class optimization if available ...
  // If this optimization is not available, the result will be wrong if T is
  // a one byte non-empty class.
  template <class T>
  class IsEmptyClass
  {
  private:
    template <class TT>
    class Test : public TT
    {
      long l;
    };

  public:
    enum
    {
      result = ((sizeof(T) == sizeof(ConstantType<T>)) && (sizeof(Test<T>) == sizeof(Test<ConstantType<T>>)))
    };
    typedef IsTrueT<result> Result;
  };

  // Type Pair : just like std::pair for types ...
  template <class A, class B>
  struct TypePair
  {
    typedef A First;
    typedef B Second;
  };

  // TypeTriplet : just like std::pair but for 3 types ...
  template <class A, class B, class C>
  struct TypeTriplet
  {
    typedef A First;
    typedef B Second;
    typedef C Third;
  };

  // findFirstHigher ... just an alias because upper_bound and especially lower_bound are confusing ;)
  template <class ForwardIterator, class T>
  ForwardIterator findFirstHigher(ForwardIterator start, ForwardIterator stop, const T &val)
  {
    return std::upper_bound(start, stop, val);
  }

  template <class ForwardIterator, class T, class Compare>
  ForwardIterator findFirstHigher(ForwardIterator start, ForwardIterator stop, const T &val,
                                  Compare comp)
  {
    return std::upper_bound(start, stop, val, comp);
  }

  template <class ForwardIterator, class T>
  ForwardIterator findFirstHigherOrEqual(ForwardIterator start, ForwardIterator stop, const T &val)
  {
    return std::lower_bound(start, stop, val);
  }

  template <class ForwardIterator, class T, class Compare>
  ForwardIterator findFirstHigherOrEqual(ForwardIterator start, ForwardIterator stop, const T &val, Compare comp)
  {
    return std::lower_bound(start, stop, val, comp);
  }

  // MINIMAL_INTEGER_TYPE (to hold (1<<N) values)
  template <int N, bool SIGNED>
  struct MinimalIntegerType : MinimalIntegerType<N - 1, SIGNED>
  {
  };

  template <>
  struct MinimalIntegerType<0, false> : ConstantType<unsigned char>
  {
  };
  template <>
  struct MinimalIntegerType<0, true> : ConstantType<char>
  {
  };
  template <>
  struct MinimalIntegerType<9, false> : ConstantType<unsigned short>
  {
  };
  template <>
  struct MinimalIntegerType<8, true> : ConstantType<short>
  {
  };
  template <>
  struct MinimalIntegerType<17, false> : ConstantType<unsigned int>
  {
  };
  template <>
  struct MinimalIntegerType<16, true> : ConstantType<int>
  {
  };

  typedef IF_<sizeof(long) == 8, long, long long>::Result
      LongInt64;
  typedef IF_<sizeof(unsigned long) == 8, unsigned long, unsigned long long>::Result
      ULongInt64;

  typedef IF_<sizeof(long long) == 16,
              long long,
              void>::Result
      LongInt128;
  typedef IF_<sizeof(unsigned long long) == 16,
              unsigned long long,
              unsigned long long>::Result
      ULongInt128;

  template <>
  struct MinimalIntegerType<33, false> : ConstantType<ULongInt64>
  {
  };
  template <>
  struct MinimalIntegerType<32, true> : ConstantType<LongInt64>
  {
  };

  // FIXME: define a 128bit integer type ??
  template <>
  struct MinimalIntegerType<65, false> : ConstantType<ULongInt128>
  {
  };
  template <>
  struct MinimalIntegerType<64, true> : ConstantType<LongInt128>
  {
  };

  // **** COUNT_BITS

  template <unsigned long long N>
  struct CountBits : ConstantValue<1 + CountBits<(N >> 1)>::value>
  {
  };
  template <>
  struct CountBits<0> : ConstantValue<0>
  {
  };
  template <>
  struct CountBits<1> : ConstantValue<1>
  {
  };

  // **** Optimal sorts

  template <class T, int N, class Reverse = IsDisabled>
  struct OptimalSortT
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      printf("Optimal sort not implemented for %d elements.\n", N);
      exit(-1);
    }
  };

  template <class T, class R>
  struct OptimalSortT<T, 8, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[0], arr[1], R());
      swap(arr[2], arr[3], R());
      swap(arr[0], arr[2], R());
      swap(arr[1], arr[3], R());
      swap(arr[1], arr[2], R());
      swap(arr[4], arr[5], R());
      swap(arr[6], arr[7], R());
      swap(arr[4], arr[6], R());
      swap(arr[5], arr[7], R());
      swap(arr[5], arr[6], R());
      swap(arr[0], arr[4], R());
      swap(arr[1], arr[5], R());
      swap(arr[1], arr[4], R());
      swap(arr[2], arr[6], R());
      swap(arr[3], arr[7], R());
      swap(arr[3], arr[6], R());
      swap(arr[2], arr[4], R());
      swap(arr[3], arr[5], R());
      swap(arr[3], arr[4], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 7, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[1], arr[2], R());
      swap(arr[0], arr[2], R());
      swap(arr[0], arr[1], R());
      swap(arr[3], arr[4], R());
      swap(arr[5], arr[6], R());
      swap(arr[3], arr[5], R());
      swap(arr[4], arr[6], R());
      swap(arr[4], arr[5], R());
      swap(arr[0], arr[4], R());
      swap(arr[0], arr[3], R());
      swap(arr[1], arr[5], R());
      swap(arr[2], arr[6], R());
      swap(arr[2], arr[5], R());
      swap(arr[1], arr[3], R());
      swap(arr[2], arr[4], R());
      swap(arr[2], arr[3], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 6, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[1], arr[2], R());
      swap(arr[0], arr[2], R());
      swap(arr[0], arr[1], R());
      swap(arr[4], arr[5], R());
      swap(arr[3], arr[5], R());
      swap(arr[3], arr[4], R());
      swap(arr[0], arr[3], R());
      swap(arr[1], arr[4], R());
      swap(arr[2], arr[5], R());
      swap(arr[2], arr[4], R());
      swap(arr[1], arr[3], R());
      swap(arr[2], arr[3], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 5, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[0], arr[1], R());
      swap(arr[3], arr[4], R());
      swap(arr[2], arr[4], R());
      swap(arr[2], arr[3], R());
      swap(arr[0], arr[3], R());
      swap(arr[0], arr[2], R());
      swap(arr[1], arr[4], R());
      swap(arr[1], arr[3], R());
      swap(arr[1], arr[2], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 4, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[0], arr[1], R());
      swap(arr[2], arr[3], R());
      swap(arr[0], arr[2], R());
      swap(arr[1], arr[3], R());
      swap(arr[1], arr[2], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 3, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[1], arr[2], R());
      swap(arr[0], arr[2], R());
      swap(arr[0], arr[1], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 2, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
      swap(arr[0], arr[1], R());
    }
  };
  template <class T, class R>
  struct OptimalSortT<T, 1, R>
  {
    static void swap(T &a, T &b, IsEnabled)
    {
      if (a < b)
        std::swap(a, b);
    }
    static void swap(T &a, T &b, IsDisabled)
    {
      if (a > b)
        std::swap(a, b);
    }
    static void sort(T *arr)
    {
    }
  };

  // math stuff
  // dot product
  template <int D, class T>
  struct DotProduct
  {
    static inline T compute(const T *a, const T *b)
    {
      return DotProduct<D - 1, T>::dot(a + 1, b + 1) + *a * *b;
    }
  };
  template <class T>
  struct DotProduct<1, T>
  {
    static inline T compute(const T *a, const T *b)
    {
      return *a * *b;
    }
  };
  template <class T>
  struct DotProduct<0, T>
  {
    static inline T compute(const T *a, const T *b)
    {
      return 0;
    }
  };

  // distance 2
  template <class T>
  static inline T square(const T a)
  {
    return a * a;
  }

  template <int D, class T>
  struct Distance2
  {
    static inline T compute(const T *a, const T *b)
    {
      return Distance2<D - 1, T>::compute(a + 1, b + 1) + square(*b - *a);
    }
  };
  template <class T>
  struct Distance2<1, T>
  {
    static inline T compute(const T *a, const T *b)
    {
      return square(*b - *a);
    }
  };
  template <class T>
  struct Distance2<0, T>
  {
    static inline T compute(const T *a, const T *b)
    {
      return 0;
    }
  };

  // norm 2
  template <int D, class T>
  struct Norm2
  {
    static inline T compute(const T *a)
    {
      return Norm2<D - 1, T>::compute(a + 1) + square(*a);
    }
  };
  template <class T>
  struct Norm2<1, T>
  {
    static inline T compute(const T *a)
    {
      return square(*a);
    }
  };
  template <class T>
  struct Norm2<0, T>
  {
    static inline T compute(const T *a)
    {
      return 0;
    }
  };

  // **** OTHER STUFF THAT MAY BE USEFULL AT SOME POINT

  // This is used to handle indices to iterate within a subarray
  // delta[i] = stride[i+1]-(dm[i]*stride[i]);
  template <int ND, typename T1, typename T2, typename T3>
  inline void getNext(T1 w[], T2 &k, const T1 dm[], const T3 delta[])
  {
    ++w[0];
    if (w[0] >= dm[0])
    {
      k += delta[0];
      if (ND == 1)
        return;
      w[0] = 0;
      ++w[1];
      if (w[1] >= dm[1])
      {
        k += delta[1];
        if (ND == 2)
          return;
        w[1] = 0;
        ++w[2];
        if (w[2] >= dm[2])
        {
          k += delta[2];
          if (ND == 3)
            return;
          w[2] = 0;
          ++w[3];
          if (w[3] >= dm[3])
          {
            k += delta[3];
            if (ND == 4)
              return;
            w[3] = 0;
            ++w[4];
            if (w[4] >= dm[4])
            {
              k += delta[4];
              if (ND == 5)
                return;
              w[4] = 0;
              ++w[5];
              if (w[5] >= dm[5])
              {
                k += delta[5];
                if (ND == 6)
                  return;
                w[5] = 0;
                ++w[6];

                int ct = 6;
                for (;;)
                {
                  if (w[ct] >= dm[ct])
                  {
                    k += delta[ct];
                    if (ND == ct + 1)
                      return;
                    w[ct] = 0;
                    ++w[++ct];
                  }
                  else
                    return;
                }
              }
            }
          }
        }
      }
    }
  }

  template <int ND, typename T>
  inline void getNext(T w[], const T dm[])
  {
    ++w[0];
    if (w[0] >= dm[0])
    {
      if (ND == 1)
        return;
      w[0] = 0;
      ++w[1];
      if (w[1] >= dm[1])
      {
        if (ND == 2)
          return;
        w[1] = 0;
        ++w[2];
        if (w[2] >= dm[2])
        {
          if (ND == 3)
            return;
          w[2] = 0;
          ++w[3];
          if (w[3] >= dm[3])
          {
            if (ND == 4)
              return;
            w[3] = 0;
            ++w[4];
            if (w[4] >= dm[4])
            {
              if (ND == 5)
                return;
              w[4] = 0;
              ++w[5];
              if (w[5] >= dm[5])
              {
                if (ND == 6)
                  return;
                w[5] = 0;
                ++w[6];

                int ct = 6;
                for (;;)
                {
                  if (w[ct] >= dm[ct])
                  {
                    if (ND == ct + 1)
                      return;
                    w[ct] = 0;
                    ++w[++ct];
                  }
                  else
                    return;
                }
              }
            }
          }
        }
      }
    }
  }
  /*
  template<>
  inline void getNext<1>(int w[], long &k, const int dm[], const long strd[])
  {
    return;
  }
  */
  template <int ND, typename L>
  inline void index2coords(L index, L coords[ND], const L stride[ND + 1])
  {
    coords[ND - 1] = index / stride[ND];
    for (int i = ND - 1; i > 0; ++i)
    {
      index -= coords[i] * stride[i + 1];
      coords[i - 1] = index / stride[i + 1];
    }
  }

  /*
  template <int A, int B>
  struct IsEqual
  {
  public:
    static const bool result=(A==B);
  };
  */

  // finds value v in C with a tolerance 'tolerance' in log(N) time
  // NOTE: C must be sorted
  template <class container>
  typename container::const_iterator findValue(const container &C, const typename container::value_type &v, double tolerance = 1.E-10)
  {
    typedef typename container::const_iterator const_iteratorT;
    const_iteratorT pos_it = std::lower_bound(C.begin(), C.end(), v);

    if ((*pos_it) != v)
    {
      if (fabs(((*pos_it) - v) / ((*pos_it) + v)) > tolerance)
      {
        pos_it--;
        if ((*pos_it) != v)
        {
          if (fabs(((*pos_it) - v) / ((*pos_it) + v)) > tolerance)
            return C.end();
        }
      }
    }
    return pos_it;
  }

  // **** Set denormal numbers to 0
  template <class T>
  void denormalsToZero(T &val, IsTrue)
  {
    if (std::fabs(val) < std::numeric_limits<T>::min())
      val = 0;
  }
  template <class T>
  void denormalsToZero(T &val, IsFalse)
  {
  }

  template <class T>
  void denormalsToZero(T &val)
  {
    typedef typename IsTrueT<SameType<T, double>::value || SameType<T, float>::value>::Result
        IsFloat;

    denormalsToZero(val, IsFloat());
  }

  template <class T>
  void denormalsToZero(T *val)
  {
    denormalsToZero(*val);
  }

};

/** \}*/
#include "../../internal/namespace.footer"
#endif
