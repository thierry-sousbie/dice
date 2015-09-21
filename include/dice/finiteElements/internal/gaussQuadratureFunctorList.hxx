#ifndef __GAUSS_QUADRATURE_FUNCTOR_LIST_HXX__
#define __GAUSS_QUADRATURE_FUNCTOR_LIST_HXX__

#include "../../internal/namespace.header"

namespace internal {
  namespace gaussQuadrature {

    class EmptyFunctor {};

    template <class F1 = EmptyFunctor, class F2 = EmptyFunctor, class F3 = EmptyFunctor,
	      class F4 = EmptyFunctor, class F5 = EmptyFunctor, class F6 = EmptyFunctor>
    class FunctorList;

    template <>
    class FunctorList<EmptyFunctor,EmptyFunctor,EmptyFunctor,
		      EmptyFunctor,EmptyFunctor,EmptyFunctor>
    {
    public:
      FunctorList(EmptyFunctor *f1=(EmptyFunctor*)(NULL), 
		  EmptyFunctor *f2=(EmptyFunctor*)(NULL), 
		  EmptyFunctor *f3=(EmptyFunctor*)(NULL), 
		  EmptyFunctor *f4=(EmptyFunctor*)(NULL), 
		  EmptyFunctor *f5=(EmptyFunctor*)(NULL), 
		  EmptyFunctor *f6=(EmptyFunctor*)(NULL))
      {}
      static const int COUNT = 0;    
    };

    template <class F1, class F2, class F3, class F4, class F5, class F6>
    class FunctorList : public FunctorList<F2,F3,F4,F5,F6>
    { 
    public:
      typedef FunctorList<F2,F3,F4,F5,F6> Next;
      typedef F1 Type;

      static const int COUNT=(Next::Count+1);

      FunctorList(F1 *f1=(F1*)(NULL), 
		  F2 *f2=(F2*)(NULL), 
		  F3 *f3=(F3*)(NULL), 
		  F4 *f4=(F4*)(NULL), 
		  F5 *f5=(F5*)(NULL), 
		  F6 *f6=(F6*)(NULL)):
	Next(f2,f3,f4,f5,f6),
	f_(f1)
      {}  
 
    
      Type *getFunctorPtr() {return f_;}

    private:
      Type *f_;
    };

    
    template <class F1=EmptyFunctor, class F2=EmptyFunctor, class F3=EmptyFunctor,
	      class F4=EmptyFunctor, class F5=EmptyFunctor, class F6=EmptyFunctor>
    FunctorList<F1,F2,F3,F4,F5,F6> 
    makeFunctorList(F1 *f1=(F1*)(NULL), F2 *f2=(F1*)(NULL), F3 *f3=(F1*)(NULL),
		    F4 *f4=(F1*)(NULL), F5 *f5=(F1*)(NULL), F6 *f6=(F1*)(NULL))
    {
      return FunctorList<F1,F2,F3,F4,F5,F6>(f1,f2,f3,f4,f5,f6);
    }
    /*
    template <class F1, class F2, class F3, class F4, class F5>
    FunctorList<F1,F2,F3,F4,F5> makeFunctorList(F1 *f1, F2 *f2, F3 *f3,F4 *f4, F5 *f5)
    {
      return FunctorList<F1,F2,F3,F4,F5>(f1,f2,f3,f4,f5);
    }

    template <class F1, class F2, class F3, class F4>
    FunctorList<F1,F2,F3,F4> makeFunctorList(F1 *f1, F2 *f2, F3 *f3, F4 *f4)
    {
      return FunctorList<F1,F2,F3,F4>(f1,f2,f3,f4);
    }

    template <class F1, class F2, class F3>
    FunctorList<F1,F2,F3> makeFunctorList(F1 *f1, F2 *f2, F3 *f3)
    {
      return FunctorList<F1,F2,F3>(f1,f2,f3);
    }

    template <class F1, class F2>
    FunctorList<F1,F2> makeFunctorList(F1 *f1, F2 *f2)
    {
      return FunctorList<F1,F2>(f1,f2);
    }

    template <class F1>
    FunctorList<F1> makeFunctorList(F1 *f1)
    {
      return FunctorList<F1>(f1);
    }
    */
    

  }
}

#include "../../internal/namespace.footer"
#endif
