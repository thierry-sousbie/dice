#ifndef __LOCAL_REGULAR_GRID_QUADRATURE_FUNCTOR_ADAPTER_HXX__
#define __LOCAL_REGULAR_GRID_QUADRATURE_FUNCTOR_ADAPTER_HXX__

#include "../../internal/namespace.header"

namespace internal {
  namespace lrg {
       
    template <class F, int ND>
    class QFunctorAdapter
    {
    public:

      static const int FE_ORDER = F::FE_ORDER;
      static const fETypeE::Type FE_TYPE = F::FE_TYPE;
      static const int QUADRATURE_DEGREE = F::QUADRATURE_DEGREE;
      
      QFunctorAdapter(F &f, int fieldIndex):
	functor(f)
      {
	fid=fieldIndex;
      }

      template <class I>
      void set(I &it)
      {
	it.coordAndSize(coords,delta);	
      }
      
      double operator()(double a, double b)
      {
	return 
	  functor(fid,coords[0]+delta[0]*a,coords[1]+delta[1]*b);	
      }

      double operator()(double a, double b, double c)
      {
	return 
	  functor(fid,coords[0]+delta[0]*a,coords[1]+delta[1]*b,coords[2]+delta[2]*c);	
      }
      
    private:
      int fid;
      F functor;      
      double delta[ND];
      double coords[ND];
    };
    /*
    template <ND>
    class QFunctorAdapter<hlp::EmptyObject>
    {};

    template <class FL, class IT>
    static void setAdaptedFunctorListHelper(FL &fl, IT &it, hlp::ConstantValue<false>)
    {}
  
    template <class FL, class IT>
    static void setAdaptedFunctorListHelper(FL &fl, IT &it, hlp::ConstantValue<true>)
    {
      //typedef typename FL::Next NextList;
      typedef typename hlp::ConstantValue< (FL::INDEX>0) > Continue;

      fl.getObject().set(it);      
      setFunctorListHelper(fl.getNext(),mesh,s,Continue());
    }
       
    template <class FL, class IT>
    void setAdaptedFunctorList(FL &fl, IT &it)
    {
      typedef typename hlp::ConstantValue< (FL::SIZE>0) > Continue;
      setFunctorListHelper(fl,mesh,s,Continue());
    }

    
    template <int ND, 
	      class F1=EmptyObject, class F2=EmptyObject, class F3=EmptyObject,
	      class F4=EmptyObject, class F5=EmptyObject, class F6=EmptyObject,
	      class F7=EmptyObject, class F8=EmptyObject, class F9=EmptyObject>
    ObjectListT<QFunctorAdapter<F1,ND>,QFunctorAdapter<F2,ND>,QFunctorAdapter<F3,ND>,
		QFunctorAdapter<F4,ND>,QFunctorAdapter<F5,ND>,QFunctorAdapter<F6,ND>,
		QFunctorAdapter<F7,ND>,QFunctorAdapter<F8,ND>,QFunctorAdapter<F9,ND> > 
    makeAdaptedObjectListHelper(int fieldIndex,
				const F1 &f1=F1(), const F2 &f2=F2(), const F3 &f3=F3(),
				const F4 &f4=F4(), const F5 &f5=F5(), const F6 &f6=F6(),
				const F7 &f7=F7(), const F8 &f8=F8(), const F9 &f9=F9())
    {
      return hlp::makeObjectList(QFunctorAdapter<F1,ND>(f1,fieldIndex),
				 QFunctorAdapter<F2,ND>(f2,fieldIndex),
				 QFunctorAdapter<F3,ND>(f3,fieldIndex),
				 QFunctorAdapter<F4,ND>(f4,fieldIndex),
				 QFunctorAdapter<F5,ND>(f5,fieldIndex),
				 QFunctorAdapter<F6,ND>(f6,fieldIndex),
				 QFunctorAdapter<F7,ND>(f7,fieldIndex),
				 QFunctorAdapter<F8,ND>(f8,fieldIndex),
				 QFunctorAdapter<F9,ND>(f9,fieldIndex));		
    }  

    template <int ND, class FL>
    ObjectListT<
      QFunctorAdapter<typename FL::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::TYPE,ND>,
      QFunctorAdapter<typename FL::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::NEXT::TYPE,ND>
      > 
    makeAdaptedObjectList(FL &fl, int fieldIndex)
    {
      return makeAdatpedObjectListHelper<ND>
	(fieldIndex,
	 fl::getObject(),
	 fl::getNext()::getObject(),
	 fl::getNext()::getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getNext()::getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getNext()::getNext()::getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getNext()::getNext()::getNext()::
	 getNext()::getObject(),
	 fl::getNext()::getNext()::getNext()::getNext()::getNext()::getNext()::
	 getNext()::getNext()::getObject());
    }
    */
  }
}
#include "../../internal/namespace.footer"
#endif
