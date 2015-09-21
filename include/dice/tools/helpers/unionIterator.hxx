#ifndef __UNION_ITERATOR_HXX__
#define __UNION_ITERATOR_HXX__

#include <iterator>


#include "../../internal/namespace.header"

template <class IA, class IB>
class UnionIterator2T
{
public:
  typedef typename IA::iterator_category iterator_category;
  typedef UnionIterator2T<IA,IB> self_type;

  typedef typename IA::container_type container_type;
  typedef typename IA::Type Type;
  typedef typename IA::value_type value_type;
  typedef typename IA::pointer pointer;
  typedef typename IA::reference reference;
  typedef typename IA::const_pointer const_pointer;
  typedef typename IA::const_reference const_reference;
  typedef typename IA::difference_type difference_type;

  UnionIterator2T(const IA &beginA_,const IA &endA_,const IB &beginB_,const IB &endB_):
    curA(beginA_),curB(beginB_),
    endA(endA_),endB(endB_)
  { 
  }

  UnionIterator2T(const IA &endA_,const IB &endB_):
    curA(endA_),curB(endB_),
    endA(endA_),endB(endB_)    
  {
  }

  virtual ~UnionIterator2T()
  {   
  }

 value_type operator->() const
  {
    if (curA!=endA)
      return curA.operator->();
    else
      return curB.operator->();
  }

  value_type operator*() const
  {    
    if (curA!=endA)
      return curA.operator*();
    else
      return curB.operator*();
  }

  self_type &operator++() // stride must be lower than the minimal size of a chunk !
  {
    if (curA!=endA)
      curA.operator++();	
    else
      curB.operator++();
      		      
    return *this;
  }

  const self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }

   bool operator==(const self_type& r) const
  {return (curA==r.curA)&&(curB==r.curB);}

  bool operator!=(const self_type& r) const
  {return (curA!=r.curA)||(curB!=r.curB);}

protected:
  IA curA;
  IB curB;
  IA endA;
  IB endB;
};


template <class IA, class IB, class IC>
class UnionIterator3T
{
public:
  typedef typename IA::iterator_category iterator_category;
  typedef UnionIterator3T<IA,IB,IC> self_type;

  typedef typename IA::container_type container_type;
  typedef typename IA::Type Type;
  typedef typename IA::value_type value_type;
  typedef typename IA::pointer pointer;
  typedef typename IA::reference reference;
  typedef typename IA::const_pointer const_pointer;
  typedef typename IA::const_reference const_reference;
  typedef typename IA::difference_type difference_type;

  UnionIterator3T(const IA &beginA_,const IA &endA_,const IB &beginB_,const IB &endB_,const IC &beginC_,const IC &endC_):
    curA(beginA_),curB(beginB_),curC(beginC_),
    endA(endA_),endB(endB_),endC(endC_)
  { 
  }

  UnionIterator3T(const IA &endA_,const IB &endB_,const IC &endC_):
    curA(endA_),curB(endB_),curC(endC_),
    endA(endA_),endB(endB_),endC(endC_)
  {
  }

  virtual ~UnionIterator3T()
  {   
  }

 value_type operator->() const
  {
    if (curA!=endA)
      return curA.operator->();
    else if (curB!=endB)
      return curB.operator->();
    else 
      return curC.operator->();    
  }

  value_type operator*() const
  {    
    if (curA!=endA)
      return curA.operator*();
    else if (curB!=endB)
      return curB.operator*();
    else 
      return curC.operator*();
  }

  self_type &operator++() // stride must be lower that the minimal size of a chunk !
  {
    if (curA!=endA)
      curA.operator++();	
    else if (curB!=endB)
      curB.operator++();
    else 
      curC.operator++();
      		      
    return *this;
  }

  const self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }

   bool operator==(const self_type& r) const
  {return (curA==r.curA)&&(curB==r.curB)&&(curC==r.curC);}

  bool operator!=(const self_type& r) const
  {return (curA!=r.curA)||(curB!=r.curB)||(curC!=r.curC);}

protected:
  IA curA;
  IB curB;
  IC curC;
  IA endA;
  IB endB;
  IC endC;
};

#include "../../internal/namespace.footer"
#endif
