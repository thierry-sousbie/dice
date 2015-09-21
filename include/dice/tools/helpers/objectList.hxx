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
    typedef ObjectListT<EmptyObject,EmptyObject,EmptyObject,
			EmptyObject,EmptyObject,EmptyObject,
			EmptyObject,EmptyObject,EmptyObject> MyType;
    typedef MyType Next;
    typedef EmptyObject Type;

    ObjectListT(const EmptyObject &f1=EmptyObject(), 
		const EmptyObject &f2=EmptyObject(), 
		const EmptyObject &f3=EmptyObject(), 
		const EmptyObject &f4=EmptyObject(), 
		const EmptyObject &f5=EmptyObject(), 
		const EmptyObject &f6=EmptyObject(), 
		const EmptyObject &f7=EmptyObject(), 
		const EmptyObject &f8=EmptyObject(), 
		const EmptyObject &f9=EmptyObject())
    {}
    static const int INDEX = -1;     
    static const int SIZE = INDEX+1;

    Type &getObject() {return obj;}
    const Type &getObject() const {return obj;}

    Next &getNext() {return static_cast<Next &>(*this);}
    const Next &getNext() const {return static_cast<const Next &>(*this);}

  private:    
    Type obj;
  };

  template <class F1,class F2,class F3,
	    class F4,class F5,class F6,
	    class F7,class F8,class F9>
  class ObjectListT : public ObjectListT<F2,F3,F4,F5,F6,F7,F8,F9>
  { 
  public:
    typedef ObjectListT<F1,F2,F3,F4,F5,F6,F7,F8,F9> MyType;
    typedef ObjectListT<F2,F3,F4,F5,F6,F7,F8,F9> Next;
    typedef F1 Type;

    static const int INDEX=(Next::INDEX+1);
    static const int SIZE=INDEX+1;
    
    ObjectListT(const F1 &f1=F1(), 
		const F2 &f2=F2(), 
		const F3 &f3=F3(), 
		const F4 &f4=F4(), 
		const F5 &f5=F5(), 
		const F6 &f6=F6(), 
		const F7 &f7=F7(), 
		const F8 &f8=F8(), 
		const F9 &f9=F9()):
      Next(f2,f3,f4,f5,f6,f7,f8,f9),
      obj(f1)
    {}  
 
    
    Type &getObject() {return obj;}
    const Type &getObject() const {return obj;}

    Next &getNext() {return static_cast<Next &>(*this);}
    const Next &getNext() const {return static_cast<const Next &>(*this);}

  private:
    Type obj;
  };

    
  template <class F1=EmptyObject, class F2=EmptyObject, class F3=EmptyObject,
	    class F4=EmptyObject, class F5=EmptyObject, class F6=EmptyObject,
	    class F7=EmptyObject, class F8=EmptyObject, class F9=EmptyObject>
  ObjectListT<F1,F2,F3,F4,F5,F6,F7,F8,F9> 
  makeObjectList(const F1 &f1=F1(), const F2 &f2=F2(), const F3 &f3=F3(),
		 const F4 &f4=F4(), const F5 &f5=F5(), const F6 &f6=F6(),
		 const F7 &f7=F7(), const F8 &f8=F8(), const F9 &f9=F9())
  {
    return ObjectListT<F1,F2,F3,F4,F5,F6,F7,F8,F9>(f1,f2,f3,f4,f5,f6,f7,f8,f9);
  }  

}

#include "../../internal/namespace.footer"
#endif
