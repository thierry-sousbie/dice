#ifndef __MY_HANDLE_HXX__
#define __MY_HANDLE_HXX__

/**
 * @file 
 * @brief  A Handle class used to unify the interface of objects and pointers.
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 *   \{
 */

/**
 * \class HandleT
 * \brief defines a class HandleT that can hold an object or a pointer to an object while 
 * giving them the same interface as a pointer. This can be usefull to hide the underlying
 * nature of the handled object.
 */

template <class C>
class HandleT
{
public:
  typedef C HandledType;

  HandleT(const C &obj):c(obj)
  {}

  HandleT(const C* obj):c(*obj)
  {}

  HandleT()
  {
    c=C();
  }

  C &operator*()
  {
    return c;
  }
  
  const C &operator*() const
  {
    return c;
  }

  C *operator->()
  {
    return &c;
  }

  const C *operator->() const
  {
    return &c;
  }

  void set(const C &obj)
  {
    c = obj;
  }

  void set(const C *obj)
  {
    c = (*obj);
  }

  HandledType *asPointer()
  {
    return &c;
  }

  const HandledType *asPointer() const
  {
    return &c;
  }

private:
  C c;  
};

template <class C>
class HandleT<C*>
{
public:
  typedef C* HandledType;

  HandleT(C *obj):c(obj)
  {}

  HandleT():c(NULL)
  {}

  C &operator*()
  {
    return *c;
  }

  const C &operator*() const
  {
    return *c;
  }

  C *operator->()
  {
    return c;
  }

  const C *operator->() const
  {
    return c;
  }

  void set(C *obj)
  {
    c = obj;
  }

  HandledType *asPointer()
  {
    return c;
  }

  const HandledType *asPointer() const
  {
    return c;
  }

private:
  C* c;
};

template <class C>
class ConstHandleT
{
public:
  typedef const C HandledType;

  ConstHandleT(const C &obj):c(obj)
  {}

  ConstHandleT(const C *obj):c(*obj)
  {}

  ConstHandleT(const HandleT<HandledType> obj):c(*obj)
  {}

  ConstHandleT():c(C())
  {}

  const C &operator*() const
  {
    return c;
  }

  const C *operator->() const
  {
    return &c;
  }

  void set(const C &obj)
  {
    c = obj;
  }

  void set(const C *obj)
  {
    c = *obj;
  }

  //! returns a pointer to the handled object
  const HandledType *asPointer() const
  {
    return &c;
  }

private:
  mutable C c;  
};

template <class C>
class ConstHandleT<C*>
{
public:
  typedef const C* HandledType;

  ConstHandleT(const C *obj):c(obj)
  {}

  ConstHandleT(const HandleT<HandledType> obj):c(obj.c)
  {}
  
  ConstHandleT():c(NULL)
  {}

  C &operator*()
  {
    return *c;
  }

  const C &operator*() const
  {
    return *c;
  }

  const C *operator->() const
  {
    return c;
  }

  void set(const C *obj)
  {
    c = obj;
  }

  //! returns a pointer to the handled object
  const HandledType *asPointer() const
  {
    return c;
  }

private:
  mutable C* c;
};

/** \}*/
#include "../../internal/namespace.footer"
#endif
