#ifndef __IMP_CLASS_WITH_SIGNATURE_ITERABLE_MEMORY_POOL_HXX__
#define __IMP_CLASS_WITH_SIGNATURE_ITERABLE_MEMORY_POOL_HXX__

#include "../../../internal/namespace.header"

template <class T, class IteratorPolicy, bool IS_POD> class IterableMemoryPoolT;

namespace internal {

  // NOTE:
  // Not 100% sure what happens on destruction of ClassWithSignature
  // FIXME: check if destructors are called correctly, whether they
  // are virtual or not ...
  // Will ~ClassWithSignature() be called twice ?

  

  // POD data are serialized as is
  // Note that big endian / little endian conversion will not work properly !
  template <class T, class IteratorPolicy, bool IS_POD = false>
  class ClassWithSignature : public T
  {
  protected:
    friend class IterableMemoryPoolT<T,IteratorPolicy,IS_POD>;
    bool iAmFree;
    void setFree() {iAmFree=true;}
    void setUsed() {iAmFree=false;}
  public:
    typedef T Base;
    typedef ClassWithSignature<T,IteratorPolicy,IS_POD> MyType;
    bool isFree() const {return iAmFree;}
    bool isUsed() const {return !iAmFree;}
    ClassWithSignature():T()
    {
      setUsed();
    }
  
    ~ClassWithSignature()
    {
      if (isUsed()) this->T::~T();
      setFree();
    }
  
    //NB: We skip reading / writing iAmFree, as any unserialized element is always used !
    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {
      writer->write(static_cast<const T*>(me));
      //char tmp=(char)me->iAmFree;
      //writer->write(&tmp);
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      new (me) MyType(); // construct me !!!!
      reader->read(static_cast<T*>(me));  
      // char tmp;
      // reader->read(&tmp);
      // me->iAmFree=(bool)tmp;
    }
  };

  // NON POD data, must have a static void selfSerialize(W*,T*);
  template <class T, class IteratorPolicy>
  class ClassWithSignature<T,IteratorPolicy,false> : public T
  {
  protected:
    friend class IterableMemoryPoolT<T,IteratorPolicy,false>;
    bool iAmFree;
    void setFree() {iAmFree=true;}
    void setUsed() {iAmFree=false;}

  public:
    typedef T Base;
    typedef ClassWithSignature<T,IteratorPolicy,false> MyType;
    bool isFree() const {return iAmFree;}
    bool isUsed() const {return !iAmFree;}

    ClassWithSignature():T()
    {
      setUsed();
    }
 
    ~ClassWithSignature()
    {
      if (isUsed()) this->T::~T();
      setFree();
    }
  

    //NB: We skip reading / writing iAmFree, as any unserialized element is always used !
    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {
      T::selfSerialize(static_cast<const T*>(me),writer);
      // char tmp=(char)me->iAmFree;
      // writer->write(&tmp);    
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      new (me) MyType(); // construct me !!!!
      T::selfUnSerialize(static_cast<T*>(me),reader);    
      // char tmp;
      // reader->read(&tmp);
      // me->iAmFree=(bool)tmp;
    } 
  };
} // internal
#include "../../../internal/namespace.footer"
#endif
