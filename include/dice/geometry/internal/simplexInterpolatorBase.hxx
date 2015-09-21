#ifndef __SIMPLEX_INTERPOLATOR_BASE_HXX__
#define __SIMPLEX_INTERPOLATOR_BASE_HXX__

#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

namespace internal {

  template <int ND>
  struct SimplexInterpolatorBaseT
  {  
    template <class T, class T2>
    static T2 interpolate(T lambda[ND+1], const T2 value[ND+1], hlp::ConstantValue<0>)
    {
      return *value;
    }
   
    template <class T, class T2>
    static T2 interpolate(T lambda[ND+1], const T2 value[ND+1], hlp::ConstantValue<1>)
    {
      T2 result=0;
      for (int i=0;i<ND+1;++i)
	result+=lambda[i]*value[i];
      return result;
    }

    template <class T, class T2, int order>
    static T2 interpolateInside(T lambda[ND+1], const T2 value[ND+1])
    {
      lambda[ND]=0;
      for (int i=0;i<ND;++i)
	{
	  if (lambda[i]<0) lambda[i]=0; 
	  if (lambda[i]>1) lambda[i]=1;
	  lambda[ND]+=lambda[i];
	}
        
      // Force the point to be inside the simplex if it is not (computation errors)
      if (lambda[ND]>1)
	{	 	  
	  double norm_inv=(1.0L)/lambda[ND];
	  for (int i=0;i<ND;++i) lambda[i]*=norm_inv;
	  lambda[ND]=0;
	}
      else lambda[ND]=1.0L-lambda[ND];

      return interpolate(lambda,value,hlp::ConstantValue<order>());
    }
  };

}

#include "../../internal/namespace.footer"
#endif

