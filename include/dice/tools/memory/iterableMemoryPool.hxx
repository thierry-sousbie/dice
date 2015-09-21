#ifndef __ITERABLE_MEMORY_POOL_HXX__
#define __ITERABLE_MEMORY_POOL_HXX__

#include <algorithm>
#include <utility>   

#include "../../dice_globals.hxx"

#include "./memoryPool.hxx"

#include "./iterableMemoryPoolPolicies.hxx"

#include "./internal/memoryPoolIterators.hxx"
#include "./internal/classWithSignature.hxx"
#include "./internal/memoryPoolValueIterator.hxx"

#include "../sort/ompPSort.hxx"

/**
 * @file 
 * @brief Defines an iterable memory pool class. IMPORTANT : The destructor of T cannot be
 *        virtual !!!!!!!
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 *   \{
 */

template <class T, 
	  class IteratorPolicy=iteratorThreadModel::Blocks,
	  bool IS_POD = false>
class IterableMemoryPoolT : 
  public MemoryPoolT< internal::ClassWithSignature<T,IteratorPolicy,IS_POD>, false >
{
protected: 
  class SortedPointerUpdaterImpl;
  friend class SortedPointerUpdaterImpl;

  // Interface to sorted pointer updater
  template <class Impl>
  class SortedPointerUpdaterInterface
  {
    template<class,class,bool> friend class IterableMemoryPoolT;
  public:
    SortedPointerUpdaterInterface(){}
    
    template <class RT>
    RT* operator()(RT* ptr) const
    {
      return (*impl)(ptr);
    }
  protected:
    SortedPointerUpdaterInterface(Impl *implementation):
      impl(implementation)
    {}
  private:
    std::shared_ptr<Impl> impl;
  };

public:  
  static std::string classHeader() {return "iterable_memory_pool";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}
  
  typedef internal::ClassWithSignature<T,IteratorPolicy,IS_POD> SignedClass;
  typedef MemoryPoolT< internal::ClassWithSignature<T,IteratorPolicy,IS_POD> , false> Base;
  typedef IterableMemoryPoolT<T,IteratorPolicy,IS_POD> MyType;

  typedef T  Type;  
  typedef T* value_type; 

  typedef SortedPointerUpdaterInterface<SortedPointerUpdaterImpl> SortedPointerUpdater;
 
  typedef internal::MemoryPoolIteratorT<MyType,value_type,IteratorPolicy> iterator;
  typedef internal::ValueIteratorT< iterator > value_iterator;
  typedef internal::MemoryPoolIteratorT<const MyType,const value_type,IteratorPolicy> const_iterator;
  typedef internal::ValueIteratorT< const_iterator > const_value_iterator;

  typedef typename Base::UnserializedPointerUpdater UnserializedPointerUpdater;

  friend class internal::MemoryPoolIteratorT<MyType,value_type,IteratorPolicy>;
  friend class internal::MemoryPoolIteratorT<const MyType,const value_type,IteratorPolicy>;
  friend class internal::ValueIteratorT< internal::MemoryPoolIteratorT<MyType,value_type,IteratorPolicy> >;
  friend class internal::ValueIteratorT< internal::MemoryPoolIteratorT<const MyType,const value_type,IteratorPolicy> >;

  IterableMemoryPoolT(std::string elName = std::string("elements"), 
		      double allocFactor_=1.0, 
		      long nMin_=10000, 
		      long granularity_ = 1):
    Base(elName,allocFactor_,nMin_,granularity_),
    nSortedElements(0)
  {
    /*
    if (sizeof(T)<=sizeof(MemoryPoolListStruct))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("IterableMemoryPool can only hold objects larger than %ld bytes\n",sizeof(MemoryPoolListStruct));	
	exit(0);
      }
    */
    if (sizeof(T)!=sizeof(SignedClass))
      glb::console->print<LOG_DEBUG>("IterableMemoryPool makes %s class grow from %ld to %ld bytes.\n",
				     Base::elementNameStr.c_str(),sizeof(T),sizeof(SignedClass));
  }

  virtual ~IterableMemoryPoolT() 
  {
    // Because this version is iterable, we can call the destructors ...
    // FIXME: behavior is inconsistent with the memoryPool, should we let the user do it ?
    if (Base::getUsedCount())
      {
	for (iterator it=begin();it!=end();++it)
	  recycle(*it);
      }
  }


private:
  // No copies please, it would not make sense ;)  
  //IterableMemoryPoolT(){}
  IterableMemoryPoolT(const MyType &init){}
  MyType &operator=(const MyType & other){}

  class GetTrueIndexFunctor
  { 
  private:
    MyType *pool;
    const std::vector<long> &elementTrueIndex;
    bool haveRecycled;

  public:

    GetTrueIndexFunctor(MyType *pool_,
			const std::vector<long> &elementTrueIndex_,
			bool haveRecycled_):
      pool(pool_),
      elementTrueIndex(elementTrueIndex_),
      haveRecycled(haveRecycled_)
    {}

    template <class TT>
    long operator()(TT* ptr) const
    {
      if (haveRecycled)
	return elementTrueIndex[pool->indexOf(static_cast<SignedClass*>(ptr))];
      else 
	return pool->indexOf(static_cast<SignedClass*>(ptr));
    }
  };

public:
  /*
  using Base::setGranularity;
  using Base::getGranularity;
  using Base::freeChunks;
  using Base::reserve;
  using Base::getUsedCount;
  using Base::getAllocatedCount;
  using Base::getRecycledCount;
  using Base::getAllocFactor;
  using Base::setAllocFactor;
  using Base::setElementName;
  using Base::getNAllocMin;
  using Base::setNAllocMin;    
  using Base::capacity;
  using Base::size;
  
  void swap(MyType &other)
  {
    Base::swap(static_cast<Base&>(other));
  }

  long indexOf(Type* ptr) const
  {
    return Base::indexOf(static_cast<SignedClass*>(ptr));
  }
  */
  
  
  void pop(Type** ptr)
  {
    SignedClass* result;
    Base::pop(&result);
    *ptr=static_cast<Type*>(result);
  }

  void recycle(Type* ptr) 
  {
    SignedClass *result=static_cast<SignedClass *>(ptr);        
    ptr->Type::~Type();
    Base::recycle_noDestructor(result);
    result->setFree();
    if (nSortedElements>0) nSortedElements--;
  }

  void recycle(Type** ptr) 
  {
    recycle(*ptr);
    *ptr=NULL;
  }
  

  // WARNING: ptr must belong to the pool or result is undefined !
  bool isFree(const Type* ptr) const
  {
    return static_cast<const SignedClass *>(ptr)->isFree();
  }

  iterator begin()
  {
    return iterator(this,(Base::getRecycledCount()>0));
  }

  iterator end()
  {    
    return iterator(this,(Base::getRecycledCount()>0),true);
  }

  iterator begin(int delta, int stride)
  {    
    return iterator(this,delta,stride,(Base::getRecycledCount()>0));
  }

  iterator end(int delta, int stride)
  {
    return iterator(this,delta,stride,(Base::getRecycledCount()>0),true);
  }

  value_iterator value_begin()
  {
    return value_iterator(this,(Base::getRecycledCount()>0));
  }

  value_iterator value_end()
  {
    return value_iterator(this,(Base::getRecycledCount()>0),true);
  }

  value_iterator value_begin(int delta, int stride)
  {
    return value_iterator(this,delta,stride,(Base::getRecycledCount()>0));
  }

  value_iterator value_end(int delta, int stride)
  {
    return value_iterator(this,delta,stride,(Base::getRecycledCount()>0),true);
  }

  const_iterator begin() const
  {
    return const_iterator(this,(Base::getRecycledCount()>0));
  }

  const_iterator end() const
  {    
    return const_iterator(this,(Base::getRecycledCount()>0),true);
  }

  const_iterator begin(int delta, int stride) const
  {
    return const_iterator(this,delta,stride,(Base::getRecycledCount()>0));
  }

  const_iterator end(int delta, int stride) const
  {
    return const_iterator(this,delta,stride,(Base::getRecycledCount()>0),true);
  }

  const_value_iterator value_begin() const
  {
    return const_value_iterator(this,(Base::getRecycledCount()>0));
  }

  const_value_iterator value_end() const
  {
    return const_value_iterator(this,(Base::getRecycledCount()>0),true);
  }

  const_value_iterator value_begin(int delta, int stride) const
  {
    return const_value_iterator(this,delta,stride,(Base::getRecycledCount()>0));
  }

  const_value_iterator value_end(int delta, int stride) const
  {
    return const_value_iterator(this,delta,stride,(Base::getRecycledCount()>0),true);
  }

  double getIteratorWastedRatio()
  {
    if (Base::getAllocatedCount()==0) return 0;
    return 1.0F - double(Base::getUsedCount())/
      (Base::getAllocatedCount() - Base::allocatedSize.back()
       +std::distance(Base::storage.back(),Base::allocEnd));
  }

  template <class W>
  void serialize(W *writer)
  {        
    writer->writeHeader(classHeader(),classVersion());
    int elSize=sizeof(T);
    writer->write(&elSize);
    Base::serialize(writer);
  }
 
  template <class R>
  UnserializedPointerUpdater unSerialize(R *reader)
  {
    if (reader == NULL)
      return Base::unSerialize(reader);
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);
    int elSize = sizeof(T);
    int fileElSize;
    reader->read(&fileElSize);
    if (elSize != fileElSize)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Cannot load memory pool from incompatible file : element size differ (current: %dB, file: %dB)\n",elSize,fileElSize);
	exit(-1);
      }
    nSortedElements=0;
    return Base::unSerialize(reader);
  }

  template <class SimplexFunctor>
  SortedPointerUpdater sort(const SimplexFunctor &f, 
			    int nThreads=glb::num_omp_threads, 
			    bool precomputeValues=true)
  {
    typedef typename std::result_of<SimplexFunctor(T*)>::type Value;

    long nElements=Base::getUsedCount();    
    std::vector<Value> value(nElements);
  
    // elementPtr is used to recover an element pointer given its true
    // index. We cannot compute that quickly on the fly as elements within pages
    // of memory are not guaranted to be contiguous (there may be holes due to
    // recycled elements)   
    std::vector<T*> elementPtr;
    elementPtr.reserve(nElements);
    
    // elementTrueIndex is used to get the true index of an element given its
    // index in the pool (returned by Base::indexOf)
    // If we have not recycled elements, then Base::indexOf returns the
    // correct true index and we don't need elementTrueIndex
    // FIXME: ACTUALLY IT DOES NOT WORK, so set haveRecycled to true for now ;)
    bool haveRecycled = true;//(Base::getRecycledCount()>0);
    std::vector<long> elementTrueIndex;
    if (haveRecycled) //elementTrueIndex.resize(nElements+Base::getRecycledCount());  
      elementTrueIndex.resize(Base::getAllocatedCount());  
    /*
    glb::console->print<LOG_STD_ALL>("Iterating elements (%s), == %ld ? (nall=%ld,nac=%ld). %ld recycled.\n",
				     Base::getElementNameStr().c_str(),
				     nElements,Base::getAllocatedCount(),
				     nElements+Base::getRecycledCount(),
				     Base::getRecycledCount());
    */
    long elementId=0;
    const iterator it_end=end();
    for (iterator it=begin(); it!=it_end;++it,++elementId)
      {
	long id=Base::indexOf(static_cast<SignedClass*>(*it));
	elementPtr.push_back(*it);
	/*
	if ((id>=elementTrueIndex.size())||(id<0))
	  glb::console->print<LOG_STD_ALL>
	    ("Invalid ID (%s): %ld>=%ld (had %ld+%ld)\n",
	     Base::getElementNameStr().c_str(),id,elementTrueIndex.size(),
	     nElements,Base::getRecycledCount());
	else
	*/
	if (haveRecycled) elementTrueIndex[id]=elementId;
      } 
    
    // Create a functor that returns the true index of an element pointer
    MyType *pool=this;
    /*
    auto getTrueIndex = 
      [&elementTrueIndex,haveRecycled,pool](T *ptr) -> long
      {
	if (haveRecycled)
	  return elementTrueIndex[pool->indexOf(static_cast<SignedClass*>(ptr))];
	else 
	  return pool->indexOf(static_cast<SignedClass*>(ptr));
      };
    */
    // For compatibility with older compilers ...
    GetTrueIndexFunctor getTrueIndex(pool,elementTrueIndex,haveRecycled);
    
    // Precompute the value of each element
#pragma omp parallel for num_threads(nThreads) //schedule(static,1)
    for (int i=0;i<nThreads;++i)
      {
	const iterator it_end = end(i,nThreads);
	for (iterator it = begin(i,nThreads); it!=it_end; ++it)
	  {	   
	    long id=getTrueIndex(static_cast<SignedClass*>(*it));	   
	    value[id] = f(*it);
	  }
      }
    
    // Initializes new indices
    std::vector<long> newIndex(nElements);
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<nElements;++i)
      newIndex[i]=i;

    // Sort elements indices in ascending order of their corresponding values   
    
    // This Lambda function declaration makes Clang++ 3.4 crash due to a bug.     
    //ompPSort(newIndex.begin(),newIndex.end(),nThreads,
    //	     [&value](long a, long b)
    //	     {return value[a]<value[b];}
    //	     );

    struct CmpValue
    {
      const std::vector<Value> &value;
      CmpValue(const std::vector<Value> &value_):value(value_){}
      bool operator()(long a, long b) const
      {return value[a]<value[b];}
    } cmpValue(value);
    ompPSort(newIndex.begin(),newIndex.end(),nThreads,cmpValue);
    
    // Single threaded version
    // std::sort(newIndex.begin(),newIndex.end(),[&value](long a, long b)
    // 	     {return value[a]<value[b];}
    // 	     );
    
    // make sure that the memory used to store the values is actually freed
    value=std::vector<Value>();    

    // Now we can reorder elements in memory according to their newIndex
    // use std::move or memcpy ? 
    // At least with memcpy we know exactly what is happening ...
    std::vector<bool> done(newIndex.size(),false);
    for (unsigned long i = 0; i < newIndex.size(); i++) 
      {
	if ((!done[i])&&(newIndex[i]!=i)) 
	  {
	    SignedClass tmp;
	    // tmp = std::move(*static_cast<SignedClass*>(elementPtr[i]));
	    // std::move(static_cast<SignedClass*>(elementPtr[i]),
	    // 	      static_cast<SignedClass*>(elementPtr[i])+1,
	    // 	      &tmp)
	    memcpy(&tmp,
	    	   static_cast<SignedClass*>(elementPtr[i]),
	    	   sizeof(SignedClass));
	    //SignedClass tmp = *curPtrArr[i];

	    unsigned long j=i;
	    while (newIndex[j] != i) 
	      {
		// *static_cast<SignedClass*>(elementPtr[j]) =
		//   std::move(*static_cast<SignedClass*>(elementPtr[newIndex[j]]));
		// std::move(static_cast<SignedClass*>(elementPtr[newIndex[j]]),
		// 	  static_cast<SignedClass*>(elementPtr[newIndex[j]]+1),
		// 	  static_cast<SignedClass*>(elementPtr[j]));
		memcpy(static_cast<SignedClass*>(elementPtr[j]),
		       static_cast<SignedClass*>(elementPtr[newIndex[j]]),
		       sizeof(SignedClass));
		//*elementPtr[j] = *elementPtr[newIndex[j]];

		j = newIndex[j];
		done[j]=true;
	      }
	    //*static_cast<SignedClass*>(elementPtr[j]) = std::move(tmp);
	    //std::move(&tmp,(&tmp)+1,static_cast<SignedClass*>(elementPtr[j]));
	    memcpy(static_cast<SignedClass*>(elementPtr[j]),
	     	   &tmp,
	     	   sizeof(SignedClass));
	    //(*elementPtr[j]) = tmp;
	  }
      }

    // newElementPtr[i] stores the new pointer to the element that used to
    // have true index i
    std::vector<T*> newElementPtr(nElements);
#pragma omp parallel for num_threads(nThreads)
    for (long i=0;i<nElements;++i)
      newElementPtr[newIndex[i]] = elementPtr[i];    
    
    nSortedElements = getIterableElementsCount();

    SortedPointerUpdaterImpl *updt=
      new SortedPointerUpdaterImpl(newElementPtr,elementTrueIndex,this,haveRecycled);

    return SortedPointerUpdater(updt);
  }

protected:
  //typedef typename Base::Chunk Chunk;
  typedef typename Base::Storage Storage;
  typedef typename Base::StorageIterator StorageIterator;
  typedef typename Base::value_type BaseValueType;

  long nSortedElements;

  void newChunk(size_t forcedSize=0)
  {
    Base::newChunk(forcedSize);
    SignedClass *curChunk = Base::storage.back();
    for (unsigned long i=0;i<Base::allocatedSize.back();i++)
      curChunk[i].setFree();   
  }

  BaseValueType getStorageBegin(const StorageIterator &it)
  {
    //return &((**it).front());
    if (Base::storage.size()==0) return NULL;
    return *it;
  }

  BaseValueType getStorageEnd(const StorageIterator &it)
  {  
    if (Base::storage.size()==0) return NULL;

    if ( (*it)==Base::storage.back() ) 
      return Base::allocEnd;
    else
      return (*it)+Base::allocatedSize[std::distance(Base::storage.begin(),it)];//+chunksize ???
  }

  long getIterableElementsCount() const
  {
    return Base::getAllocatedCount()-Base::getNSpare();
  }

  long getSortedIterableElementsCount() const
  {
    return nSortedElements;
  }

  Storage &getStorage()
  {
    return Base::storage;
  }

  long getStorageSize(long i)
  {
    return Base::allocatedSize[i];
  }

};

template <class T, class IteratorPolicy, bool IS_POD>
class IterableMemoryPoolT<T,IteratorPolicy,IS_POD>::SortedPointerUpdaterImpl
{
  public:
  typedef T Data;
  typedef T* DataPtr;
  typedef IteratorPolicy DataIteratorPolicy;
  static const bool PODNESS = IS_POD;
  
  typedef internal::ClassWithSignature<T,IteratorPolicy,IS_POD> SignedClass;
  typedef IterableMemoryPoolT<T,IteratorPolicy,IS_POD> Pool;
  friend class IterableMemoryPoolT<T,IteratorPolicy,IS_POD>;
            
  template <class RT>
  RT* operator()(RT* ptr) const
  {
    return newElementPtr[getTrueIndex(ptr)];
  }

  SortedPointerUpdaterImpl()
  {}

private:
  // Be carefull, this constructor stills the vectors !
  SortedPointerUpdaterImpl(std::vector<T*> &newElementPtr_,
			   std::vector<long> &elementTrueIndex_,
			   Pool* pool_, bool haveRecycledElements_)
  {
    newElementPtr.swap(newElementPtr_);
    elementTrueIndex.swap(elementTrueIndex_);
    pool=pool_;
    haveRecycledElements=haveRecycledElements_;//(pool->getRecycledCount()>0);
  }

  long getTrueIndex(T *ptr) const
  {
    if (haveRecycledElements)
      return elementTrueIndex[pool->indexOf(static_cast<SignedClass*>(ptr))];
    else 
      return pool->indexOf(static_cast<SignedClass*>(ptr));
  }

  std::vector<T*> newElementPtr;
  std::vector<long> elementTrueIndex;
  Pool* pool;
  bool haveRecycledElements;
};

/** \}*/
#include "../../internal/namespace.footer"
#endif
