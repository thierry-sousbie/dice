#ifndef __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX_HXX__

#include "../../internal/namespace.header"

namespace internal {
  namespace gaussQuadrature {

    template <int NVAL> struct Permutation;

    template <>
    struct Permutation<2>
    {
      template <class F> 
      static double get(F &f, double W, double A)
      {
	return W*f(A,A);
      }

      template <class F> 
      static double get(F &f, double W, double A, double B)
      {
	return W*(f(A,B)+f(B,A));
      }
    };

    template <>
    struct Permutation<3>
    {
      template <class F> 
      static double get(F &f, double W, double A)
      {
	return W*f(A,A,A);
      }

      template <class F> 
      static double get(F &f, double W, double A, double B)
      {
	return W*(f(A,A,B)+f(A,B,A)+f(B,A,A));
      }

      template <class F> 
      static double get(F &f, double W, double A, double B, double C)
      {
	return 
	  W*(f(A,B,C)+f(C,A,B)+f(B,C,A)+
	     f(B,A,C)+f(C,B,A)+f(A,C,B));
      }
    };   

    template <>
    struct Permutation<4>
    {
      template <class F> 
      static double get(F &f, double W, double A)
      {
	return W*f(A,A,A,A);
      }

      template <class F> 
      static double get31(F &f, double W, double A, double B)
      {
	return W*(f(A,A,A,B)+f(A,A,B,A)+f(A,B,A,A)+f(B,A,A,A));
      }

      template <class F> 
      static double get22(F &f, double W, double A, double B)
      {
	return W*(f(A,A,B,B)+f(B,A,A,B)+f(B,B,A,A)+f(A,B,B,A)+f(A,B,A,B)+f(B,A,B,A));
      }

      template <class F> 
      static double get211(F &f, double W, double A, double B, double C)
      {
	return 
	  W*(f(A,A,B,C)+f(A,B,A,C)+f(A,B,C,A)+f(B,A,A,C)+f(B,A,C,A)+f(B,C,A,A)+
	     f(A,A,C,B)+f(A,C,A,B)+f(A,C,B,A)+f(C,A,A,B)+f(C,A,B,A)+f(C,B,A,A));
      }

      template <class F> 
      static double get(F &f, double W, double A, double B, double C, double D)
      {
	return 
	  W*(f(A,B,C,D) + f(A,B,D,C) +	     
	     f(A,C,B,D) + f(A,D,B,C) +
	     f(A,C,D,B) + f(A,D,C,B) +

	     f(B,A,C,D) + f(B,A,D,C) +	     
	     f(C,A,B,D) + f(D,A,B,C) +
	     f(C,A,D,B) + f(D,A,C,B) +

	     f(B,C,A,D) + f(B,D,A,C) +	     
	     f(C,B,A,D) + f(D,B,A,C) +
	     f(C,D,A,B) + f(D,C,A,B) +

	     f(B,C,D,A) + f(B,D,C,A) +	     
	     f(C,B,D,A) + f(D,B,C,A) +
	     f(C,D,B,A) + f(D,C,B,A));
      }
    };

  }
}

#include "../../internal/namespace.footer"

#include "specializedGQ_simplex2D.hxx"
#include "specializedGQ_simplex3D.hxx"

#endif
