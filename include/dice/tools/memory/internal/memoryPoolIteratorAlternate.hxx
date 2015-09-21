#ifndef __MEMORY_POOL_ITERATOR_ALTERNATE_HXX__
#define __MEMORY_POOL_ITERATOR_ALTERNATE_HXX__

#include <iterator>
#include "memoryPoolIterators.hxx"

#include "../../../internal/namespace.header"

// This iterator gives each thread one item every nThreads (i.e. the ith thread among a 
// total of  nThreads is given items with indices i+k*nThreads)

namespace internal {

  template <class C, class V>
  class MemoryPoolIteratorT<C,V,iteratorThreadModel::Alternate>
  {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef iteratorThreadModel::Alternate ThreadTraits;
    typedef MemoryPoolIteratorT<C,V,ThreadTraits> self_type;

    typedef C container_type;
    typedef typename container_type::Type Type;
    typedef V value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long difference_type;
  
    MemoryPoolIteratorT(C *c,bool checkFreed_,bool end=false):
      stride(1),checkFreed(checkFreed_),
      container(c)
    {
      stIt = container->getStorage().begin();     
      endPtr = container->getStorageEnd(stIt);//(&((**stIt).back()) + 1);

      if ((end)||(container->getUsedCount()<1)) curPtr=NULL;
      else 
	{
	  curPtr = container->getStorageBegin(stIt);//&((**stIt).front());
	  //if (((SignedClass*)curPtr)->isFree()) this->operator++();
	  if (curPtr->isFree()) this->operator++();
	}
    }

    MemoryPoolIteratorT(C *c,long delta, long stride_,bool checkFreed_, bool end=false):
      stride(stride_),checkFreed(checkFreed_),
      container(c)
    {
      stIt = container->getStorage().begin();  
      endPtr = container->getStorageEnd(stIt);//(&((**stIt).back()) + 1);
   
      if ((end)||(container->getUsedCount()<1+delta))
	curPtr = NULL;
      else 
	{
	  curPtr = container->getStorageBegin(stIt)+delta;//(&((**stIt).front()) + delta);
	  //if (((SignedClass*)curPtr)->isFree()) this->operator++();
	  if (curPtr->isFree()) this->operator++();
	}
    }

    virtual ~MemoryPoolIteratorT()
    {
    
    }

    value_type operator->() const
    {
      return static_cast<value_type>(curPtr);
    }

    inline value_type operator*() const
    {    
      return static_cast<value_type>(curPtr);
    }

    // stride must be lower than the minimal size of a chunk !
    inline self_type &operator++() 
    {
      //printf("curPtr=%ld / end=%ld\n",(long)curPtr,(long)endPtr);
      curPtr += stride;

      // FIXME::remove the comment !!!!!
      if (curPtr>=endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();

      /*
	{
	  stIt++;
	  if (stIt!=container->getStorage().end())
	    {
	      long delta=(long)(curPtr-endPtr);
	      curPtr = container->getStorageBegin(stIt)+delta;//&((**stIt).front())+delta;
	      endPtr = container->getStorageEnd(stIt);//&((**stIt).back()) + 1;
	      //if (curPtr->isFree()) this->operator++();
	    }
	  else 
	    {
	      curPtr=NULL;
	      //return *this;
	    }
	}
      */
    
      //FIXME!
      /*
      if (curPtr->isFree()) 
	return this->operator++();
      */
      return *this;
    }

    const self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    bool operator==(const self_type& r) const
    {return (curPtr==r.curPtr);}

    bool operator!=(const self_type& r) const
    {return (curPtr!=r.curPtr);}
 
  private:
    void increment()
    {
      curPtr += stride;

      if (curPtr>=endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();

      //return *this;
    }


    //self_type& operator=( const self_type& other ){}
    void advanceNewPage()
    {
      ++stIt;
      if (stIt!=container->getStorage().end())
	{
	  long delta=(long)(curPtr-endPtr);
	  curPtr = container->getStorageBegin(stIt)+delta;//&((**stIt).front())+delta;
	  endPtr = container->getStorageEnd(stIt);//&((**stIt).back()) + 1;
	  if (checkFreed)
	    if (curPtr->isFree()) increment();
	}
      else 
	{
	  curPtr=NULL;
	  //return *this;
	}
    }
    
  protected:
    typedef typename C::StorageIterator StorageIterator;
    //typedef typename C::SignedClass SignedClass;
    typedef typename C::BaseValueType BaseValueType;
    
    BaseValueType curPtr;
    BaseValueType endPtr;  
    long stride;
    long checkFreed;
    container_type *container;  
    StorageIterator stIt;
  };

} // internal

#include "../../../internal/namespace.footer"
#endif
