#ifndef __SPECIALIZED_GAUSS_QUADRATURE_CUBOID_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_CUBOID_HXX__

#include "../../internal/namespace.header"


namespace internal {
  namespace gaussQuadrature {
  
    template <int order, typename T=double>
    struct CuboidPoints;

    template <typename T>
    struct CuboidPoints<0,T>
    {
      static const int NP  = 1;
   
      static T coef(int which)
      {	
	T C[NP]={0};
	return  C[which];
      }

      static T w(int which)
      {
	return 1.0;
      }

    };
    
    template <typename T>
    struct CuboidPoints<1,T>
    {
      static const int NP  = 1;
    
      static T coef(int which)
      {	
	T C[NP]={0.5};
	
	return  C[which];
      }

      static T w(int which)
      {
	return 1.0;
      }
    };

    template <typename T>
    struct CuboidPoints<3,T>
    {
      static const int NP  = 2;
    
      static T coef(int which)
      {
	static const T A=1.0/sqrt(3.0);
	T C[NP]={0.5*(1.0-A),0.5*(1.0+A)};
	return  C[which];
      }

      static T w(int which)
      {
	return 0.5;
      }
    };
    template <typename T>
    struct CuboidPoints<2,T> :
      public CuboidPoints<3,T>
    {};

    template <typename T>
    struct CuboidPoints<5,T>
    {
      static const int NP  = 3;
    
      static T coef(int which)
      {
	static const T A=sqrt(3.0/5.0);
	T C[NP]={0.5*(1.0-A),0.5*(1.0+A),0.5};
	return  C[which];
      }

      static T w(int which)
      {
	static const T W1=5.0/18.0;
	static const T W2=8.0/18.0;
	static const T W[NP]={W1,W1,W2};
	return W[which];
      }
    };
    template <typename T>
    struct CuboidPoints<4,T> :
      public CuboidPoints<5,T>
    {};

    template <typename T>
    struct CuboidPoints<7,T>
    {
      static const int NP  = 4;
      
      static T coef(int which)
      {
	static const T A=sqrt(3.0+2.0/7.0*sqrt(6.0/5.0));
	static const T B=sqrt(3.0-2.0/7.0*sqrt(6.0/5.0));
	
	T C[NP]={0.5*(1.0-A),0.5*(1.0+A),0.5*(1.0-B),0.5*(1.0+B)};
	return  C[which];
      }

      static T w(int which)
      {
	static const T W1=0.25-1.0/12.0*sqrt(5.0/6.0);
	static const T W2=0.25+1.0/12.0*sqrt(5.0/6.0);
	static const T W[NP]={W1,W1,W2,W2};
	return W[which];
      }
    };
    template <typename T>
    struct CuboidPoints<6,T> :
      public CuboidPoints<7,T>
    {};

    template <typename T>
    struct CuboidPoints<9,T>
    {
      static const int NP  = 5;

      static T coef(int which)
      {
	static const T A=sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
	static const T B=sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
	
	T C[NP]={0.5*(1.0-A),0.5*(1.0+A),0.5*(1.0-B),0.5*(1.0+B),0.5};
	return  C[which];
      }

      static T w(int which)
      {
	static const T W1=(322.0-13.0*sqrt(70.0))/1800.0;
	static const T W2=(322.0+13.0*sqrt(70.0))/1800.0;
	static const T W3=512.0/1800.0;
	static const T W[NP]={W1,W1,W2,W2,W3};
	return W[which];
      }
    };
    template <typename T>
    struct CuboidPoints<8,T> :
      public CuboidPoints<9,T>
    {};

  }
}

#include "../../internal/namespace.footer"

#include "specializedGQ_cuboidND.hxx"

#endif
