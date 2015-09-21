#ifndef __SPECIALIZED_GAUSS_QUADRATURE_CUBOIDND_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_CUBOIDND_HXX__

#include "specializedGQ_cuboid1D.hxx"
#include "../../tools/helpers/helpers.hxx"
#include "../../internal/namespace.header"

namespace internal {
  namespace gaussQuadrature {

    // ND gauss points coordinates are just implemented as a tensor product of the 1D case
    template <class F>
    struct CuboidFunctorAdapter1
    {
      CuboidFunctorAdapter1(F& f):
	functor(f)
      {}

      void setCoef1(double u)
      {
	c1=u;
      }

      double operator()(double a)
      {
	return functor(c1,a);
      }

      F &functor;
      double c1;
    };

    template <class F>
    struct CuboidFunctorAdapter2
    {
      CuboidFunctorAdapter2(F& f):
	functor(f)
      {}
      
      void setCoef1(double u)
      {
	c1=u;
      }

      void setCoef2(double u)
      {
	c2=u;
      }

      double operator()(double a)
      {
	return functor(c1,c2,a);
      }

      F &functor;
      double c1;
      double c2;
    };

    template <int DEG>
    struct SpecializedGQT<2,DEG,fETypeE::cuboid>
    {
      typedef SpecializedGQT<1,DEG,fETypeE::cuboid> Cuboid1D;
      typedef CuboidPoints<DEG> CP;

      static const int NDIM = 2;
      static const int NP_1D  = Cuboid1D::NP;
      static const int NP  = hlp::IntPower<NP_1D,NDIM>::value;

      template <class F>
      static double get(F& functor)
      {
	double result=0;
	CuboidFunctorAdapter1<F> f(functor);

	for (int i=0;i<NP_1D;++i)
	  {	   
	    f.setCoef1(CP::coef(i));
	    result+=CP::w(i)*Cuboid1D::get(f);
	  }

	return result;
      }   
    };

    template <int DEG>
    struct SpecializedGQT<3,DEG,fETypeE::cuboid>
    {
      typedef SpecializedGQT<1,DEG,fETypeE::cuboid> Cuboid1D;
      typedef CuboidPoints<DEG> CP;

      static const int NDIM = 3;
      static const int NP_1D  = Cuboid1D::NP;
      static const int NP  = hlp::IntPower<NP_1D,NDIM>::value;

      template <class F>
      static double get(F& functor)
      {
	double result=0;
	CuboidFunctorAdapter2<F> f(functor);

	for (int i=0;i<NP_1D;++i)
	  {	   
	    double w=CP::w(i);
	    f.setCoef1(CP::coef(i));
	    for (int j=0;j<NP_1D;++j)
	      {
		f.setCoef2(CP::coef(j));
		result+=w*CP::w(j)*Cuboid1D::get(f);
	      }
	  }

	return result;
      }   
    };

  }
}


#include "../../internal/namespace.footer"
#endif
