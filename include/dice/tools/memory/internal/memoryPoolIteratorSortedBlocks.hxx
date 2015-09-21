#ifndef __MEMORY_POOL_ITERATOR_SORTED_BLOCKS_HXX__
#define __MEMORY_POOL_ITERATOR_SORTED_BLOCKS_HXX__

#include <iterator>
#include "memoryPoolIterators.hxx"

#include "../../../internal/namespace.header"

// Each threads is given two blocks (each of them contiguous): a sorted one and a non 
// sorted one. This is supposed to be used when sorting the items improves the cache 
// access pattern significantly, having two blocks (a fast sorted one and a slow 
// random one) allows for a much better load balancing than iteratorThreadModel::Blocks

namespace internal {

  template <class C, class V>
  class MemoryPoolIteratorT<C,V,iteratorThreadModel::SortedBlocks>
  {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef iteratorThreadModel::SortedBlocks ThreadTraits;
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
      /*
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
      */
      m_delta=0;
      m_stride=1;
      init(0,end);
    }
    
    MemoryPoolIteratorT(C *c,int delta, int stride,bool checkFreed_, bool end=false):
      container(c),checkFreed(checkFreed_)
    {
      m_delta=delta;
      m_stride=stride;
      init(0,end);
      /*
      long nElements = container->getSortedIterableElementsCount();      
      init(nElements,0,end);
      */

      /*
      long startShift = container->getSortedIterableElementsCount();
      long nElements = container->getIterableElementsCount()-startShift;
      init(nElements,startShift);
      */

      /*
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
      */
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
    
      if (curPtr==endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();

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
    int getCurBlockIndex() const
    {
      return (m_stride==0);
      // if (m_stride==0) return 1;
      // return 0;
    }

    void init(int blockIndex, bool end=false)
    {
      long startShift;
      long nElements;  
      
      if (m_stride==1)
	{
	  // We do not need two blocks in that case
	  //printf("Hi there ;)\n");
	  startShift=0;
	  nElements = container->getIterableElementsCount();  
	}
      else if (blockIndex==0)
	{
	  startShift=0;
	  nElements = container->getSortedIterableElementsCount();  
	}
      else
	{
	  startShift = container->getSortedIterableElementsCount();
	  nElements = container->getIterableElementsCount()-startShift;
	}
      // printf("blockId=%d, startShift=%ld, nElements=%ld\n",
      // 	     blockIndex, startShift, nElements);

      /*
      startShift=0;
      nElements = container->getIterableElementsCount();  
      */
      StorageIterator storageBegin = container->getStorage().begin();
      long count = (nElements/m_stride);
      long excess = nElements - (count*m_stride);            
      long startIndex = (count*m_delta) + ((m_delta<excess)?m_delta:excess) + startShift;
      long stopIndex = startIndex + count + ((m_delta<excess)?1:0);

      // glb::console->print<LOG_STD_ALL>("(%d/%d): startIndex=%ld stopIndex=%ld\n",
      // 				       m_delta,m_stride,startIndex,stopIndex);
      // printf("startIndex=%ld stopIndex=%ld\n",
      // 	     startIndex,stopIndex);

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

	  if ((blockIndex>0)||(m_stride==1)) m_stride=0; // This marks the second block

	  if ((end)||(container->getUsedCount()<1+m_delta))
	    curPtr = NULL;
	  else if (curPtr->isFree()) 
	    this->operator++();
	}
      else if (blockIndex==0)
	init(1,end);
      else curPtr=NULL;

      if ((blockIndex>0)||(m_stride==1)) m_stride=0; // This marks the second block

      // printf("=> curPtr = %ld\n",(long)curPtr);
      //m_stride=0;
    }


    void increment()
    {
      ++curPtr;        
      if (curPtr==endPtr) advanceNewPage();
      else if (checkFreed)
	if (curPtr->isFree()) increment();
    }
   
    void advanceNewPage()
    {
      if (endPtr==globalEndPtr)
	{
	  if (getCurBlockIndex()==0) init(1);
	  else curPtr=NULL;	    
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
	    if (curPtr->isFree()) increment();
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
    int m_delta;
    int m_stride;
    
  };

} // internal

#include "../../../internal/namespace.footer"
#endif
