#ifndef __MEMORY_POOL_VALUE_ITERATOR_HXX__
#define __MEMORY_POOL_VALUE_ITERATOR_HXX__

#include <iterator>
#include "memoryPoolIterators.hxx"

#include "../../../internal/namespace.header"

namespace internal {
  template <class T> class ValueIteratorT;

  //template <>
  template <class C, class V, class T>
  class ValueIteratorT<  MemoryPoolIteratorT<C,V,T> >
    :public MemoryPoolIteratorT<C,V,T>
  {
    typedef MemoryPoolIteratorT<C,V,T> Base;
  public:
  
    typedef ValueIteratorT<  Base > self_type;
    typedef typename Base::Type value_type;
    typedef typename Base::value_type base_value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;

    ValueIteratorT(C *c,bool end=false):
      Base(c,end)
    {

    }

    ValueIteratorT(C *c,long delta, long stride_, bool end=false):
      Base(c,delta,stride_,end)
    {

    }

    virtual ~ValueIteratorT()
    {
    
    }

    reference operator->() const
    {
      return *static_cast<base_value_type>(Base::curPtr);
    }

    reference operator*() const
    {    
      return *static_cast<base_value_type>(Base::curPtr);
    }

    self_type &operator++()
    {
      Base::operator++();
      return *this;
    }

    const self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

  };

}

#include "../../../internal/namespace.footer"
#endif
