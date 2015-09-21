#ifndef __MEMORY_POOL_ITERATOR_BLOCKS_HXX__
#define __MEMORY_POOL_ITERATOR_BLOCKS_HXX__

#include <iterator>
#include "memoryPoolIterators.hxx"

#include "../../../internal/namespace.header"

// This iterator gives each thread a contiguous block of memory.

namespace internal {

  template <class C, class V>
  class MemoryPoolIteratorT<C,V,iteratorThreadModel::Blocks>
  {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef iteratorThreadModel::Blocks ThreadTraits;
    typedef MemoryPoolIteratorT<C,V,ThreadTraits> self_type;

    typedef C container_type;
    typedef typename container_type::Type Type;
    typedef V value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long difference_type;
  
    MemoryPoolIteratorT(C *c,bool checkFreed_, bool end=false):
      container(c),checkFreed(checkFreed_)
    {
      stIt = container->getStorage().begin();     
      stItLast = stIt+container->getStorage().size()-1;
      endPtr = container->getStorageEnd(stIt);     
      globalEndPtr = container->getStorageEnd(stItLast);    

      if ((end)||(container->getUsedCount()<1)) 
	curPtr=NULL;
      else 
	{
	  curPtr = container->getStorageBegin(stIt);
	  if (curPtr->isFree()) this->operator++();
	}
    }
    
    MemoryPoolIteratorT(C *c,long delta, long stride,bool checkFreed_, bool end=false):
      container(c),checkFreed(checkFreed_)
    {
      StorageIterator storageBegin = container->getStorage().begin();
      long nElements = container->getIterableElementsCount();      
      long count = (nElements/stride);
      long excess = nElements - (count*stride);            
      long startIndex = (count*delta) + ((delta<excess)?delta:excess);
      long stopIndex = startIndex + count + ((delta<excess)?1:0);

      if (startIndex < stopIndex)
      	{
	  long startChunkIndex = 0;
	  while (startIndex >= container->getStorageSize(startChunkIndex))
	    {
	      startIndex -= container->getStorageSize(startChunkIndex);
	      stopIndex -= container->getStorageSize(startChunkIndex);
	      startChunkIndex++;
	    }
	  
	  long stopChunkIndex = startChunkIndex;
	  while (stopIndex > container->getStorageSize(stopChunkIndex))
	    {
	      stopIndex -= container->getStorageSize(stopChunkIndex);
	      stopChunkIndex++;
	    }
	  //printf("StartChunk : %ld / %ld\n",startChunkIndex,stopChunkIndex);	
      
	  stIt = storageBegin + startChunkIndex;
	  stItLast = storageBegin + stopChunkIndex;
	  curPtr = container->getStorageBegin(stIt)+startIndex;
	  globalEndPtr = container->getStorageBegin(storageBegin+stopChunkIndex)+stopIndex;
	  
	  if (stIt == stItLast)
	    endPtr = globalEndPtr;
	  else
	    endPtr = container->getStorageEnd(stIt);   

	  if ((end)||(container->getUsedCount()<1+delta))
	    curPtr = NULL;
	  else if (curPtr->isFree()) this->operator++();
	}
      else curPtr=NULL;

      //printf("Going from %ld -> %ld(%ld)\n",curPtr,endPtr,globalEndPtr);      	
    }

    virtual ~MemoryPoolIteratorT()
    {
    
    }

    value_type operator->() const
    {
      return static_cast<value_type>(curPtr);
    }

    value_type operator*() const
    {    
      return static_cast<value_type>(curPtr);
    }

    inline self_type &operator++() 
    {
      ++curPtr;  
      // FIXME::remove the comment !!!!!
      if (curPtr==endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();//this->operator++();

      /*
	{
	  if (endPtr==globalEndPtr)
	    {
	      curPtr=NULL;
	      return *this;
	    }
	  else
	    {
	      stIt++;
	      curPtr = container->getStorageBegin(stIt);

	      if (stIt==stItLast)
		endPtr = globalEndPtr;
	      else
		endPtr = container->getStorageEnd(stIt);

	      if (curPtr->isFree()) return this->operator++();
	    }
	}
      */

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
      ++curPtr;        
      if (curPtr==endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();//this->operator++();
    }
    //self_type& operator=( const self_type& other ){}
    void advanceNewPage()
    {
      if (endPtr==globalEndPtr)
	{
	  curPtr=NULL;
	  //return *this;
	}
      else
	{
	  ++stIt;
	  curPtr = container->getStorageBegin(stIt);

	  if (stIt==stItLast)
	    endPtr = globalEndPtr;
	  else
	    endPtr = container->getStorageEnd(stIt);

	  if (checkFreed)
	    if (curPtr->isFree()) increment();//this->operator++();
	}      
    }
  
 
  protected:
    typedef typename C::StorageIterator StorageIterator;
    //typedef typename C::SignedClass SignedClass;
    typedef typename C::BaseValueType BaseValueType;

    container_type *container;  
    BaseValueType curPtr;
    BaseValueType endPtr;  
    long checkFreed;
    BaseValueType globalEndPtr;
    StorageIterator stIt;
    StorageIterator stItLast;
    
  };

} // internal

#include "../../../internal/namespace.footer"
#endif
