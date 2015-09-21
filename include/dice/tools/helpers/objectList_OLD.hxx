#ifndef __HLP_OBJECT_LIST_HXX__
#define __HLP_OBJECT_LIST_HXX__

#include "../../internal/namespace.header"

namespace hlp {

  class EmptyObject {};

  template <class F1 = EmptyObject, class F2 = EmptyObject, class F3 = EmptyObject,
	    class F4 = EmptyObject, class F5 = EmptyObject, class F6 = EmptyObject,
	    class F7 = EmptyObject, class F8 = EmptyObject, class F9 = EmptyObject>
  class ObjectListT;

  template <>
  class ObjectListT<EmptyObject,EmptyObject,EmptyObject,
		    EmptyObject,EmptyObject,EmptyObject,
		    EmptyObject,EmptyObject,EmptyObject>
  {
  public:
    ObjectListT(EmptyObject *f1=(EmptyObject*)(NULL), 
		EmptyObject *f2=(EmptyObject*)(NULL), 
		EmptyObject *f3=(EmptyObject*)(NULL), 
		EmptyObject *f4=(EmptyObject*)(NULL), 
		EmptyObject *f5=(EmptyObject*)(NULL), 
		EmptyObject *f6=(EmptyObject*)(NULL), 
		EmptyObject *f7=(EmptyObject*)(NULL), 
		EmptyObject *f8=(EmptyObject*)(NULL), 
		EmptyObject *f9=(EmptyObject*)(NULL))
    {}
    static const int COUNT = 0;  

    void *getObjectPtr() {return NULL;}
  };

  template <class F1,class F2,class F3,
	    class F4,class F5,class F6,
	    class F7,class F8,class F9>
  class ObjectListT : public ObjectListT<F2,F3,F4,F5,F6,F7,F8,F9>
  { 
  public:
    typedef ObjectListT<F2,F3,F4,F5,F6,F7,F8,F9> Next;
    typedef F1 Type;

    static const int COUNT=(Next::COUNT+1);
    
    ObjectListT(F1 *f1=(F1*)(NULL), 
		F2 *f2=(F2*)(NULL), 
		F3 *f3=(F3*)(NULL), 
		F4 *f4=(F4*)(NULL), 
		F5 *f5=(F5*)(NULL), 
		F6 *f6=(F6*)(NULL), 
		F7 *f7=(F7*)(NULL), 
		F8 *f8=(F8*)(NULL), 
		F9 *f9=(F9*)(NULL)):
      Next(f2,f3,f4,f5,f6,f7,f8,f9),
      f_(f1)
    {}  
 
    
    Type *getObjectPtr() {return f_;}

  private:
    Type *f_;
  };

    
  template <class F1=EmptyObject, class F2=EmptyObject, class F3=EmptyObject,
	    class F4=EmptyObject, class F5=EmptyObject, class F6=EmptyObject,
	    class F7=EmptyObject, class F8=EmptyObject, class F9=EmptyObject>
  ObjectListT<F1,F2,F3,F4,F5,F6,F7,F8,F9> 
  makeObjectList(F1 *f1=(F1*)(NULL), F2 *f2=(F2*)(NULL), F3 *f3=(F3*)(NULL),
		 F4 *f4=(F4*)(NULL), F5 *f5=(F5*)(NULL), F6 *f6=(F6*)(NULL),
		 F7 *f7=(F7*)(NULL), F8 *f8=(F8*)(NULL), F9 *f9=(F9*)(NULL))
  {
    return ObjectListT<F1,F2,F3,F4,F5,F6,F7,F8,F9>(f1,f2,f3,f4,f5,f6,f7,f8,f9);
  }  

}

#include "../../internal/namespace.footer"
#endif
