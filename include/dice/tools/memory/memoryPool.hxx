#ifndef __MEMORY_POOL_HXX__
#define __MEMORY_POOL_HXX__

#include <stdio.h>

#include <vector>
#include <limits>
#include <memory>
#include <functional> 
#include <algorithm>
#include <iterator>

#include "../../dice_config.hxx"
#include "../../tools/IO/myIO.hxx"
#include "../../tools/helpers/helpers.hxx"
#include "../../dice_globals.hxx"

/**
 * @file 
 * @brief Defines a (non iterable) memory pool class
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

// NB: It is the users responsibility to recycle any allocated objects before 
//     destroying a pool. Failure to do so will result in the destructors not
//     being called.

/** \addtogroup TOOLS
 *   \{
 */

struct MemoryPoolListStruct 
{
  struct MemoryPoolListStruct *next;
};

template <class T, bool IS_POD = false>
class MemoryPoolT
{
protected:
  class UnserializedPointerUpdate;
  friend class UnserializedPointerUpdate;

public:   
  typedef MemoryPoolT<T,IS_POD> MyType;
  typedef T* value_type;  
  typedef T Type;
  
  typedef std::unique_ptr<UnserializedPointerUpdate> UnserializedPointerUpdater;

  // We use 64 as it is most likely the size of a cache line + we can use avx instructions on anything 
  // starting at the begining of objects (32 would be enough for avx ...)
  static const int ALIGNMENT = 64; 

  static std::string classHeader() {return "memory_pool";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  MemoryPoolT(std::string elName = std::string("elements"), 	      
	      double allocFactor_=1.0, 
	      long nMin_=10000, 
	      long granularity_ = 1): 
    elementNameStr(elName),
    allocEnd(NULL),
    nUsed(0),
    nAllocated(0),
    //nRecycled(0),
    freeData(NULL)
  {
    if (sizeof(T)<sizeof(MemoryPoolListStruct))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("MemoryPool can only hold objects larger than %ld bytes\n",sizeof(MemoryPoolListStruct)-1);	
	exit(0);
      }
    
    setAllocFactor(allocFactor_);
    setNAllocMin(nMin_);
    reserve(0); 
    setGranularity(granularity_);
  }
  
private:
  // No copies please, it would not make sense ;)  
  //MemoryPoolT(){}
  MemoryPoolT(const MyType &other);
  MyType &operator=(const MyType & other);
  
public:

  // FIXME ? : if (nUsed!=0) then the destructors of the allocated objects
  //           will NOT be called ! This is okay if the memory pool lives
  //           as long as the application, or if data does not need to be
  //           destroyed ... 
  //           We could fix that but destructor would become a O(N.log(N_free)) 
  // N.B. : it is the user's responsability to recycle allocated objects before 
  // destroying a pool !  
  virtual ~MemoryPoolT() 
  {   
    freeChunks();
  }

  void setGranularity(unsigned long g)
  {    
    granularity=(g<1)?1:g;
  }

  long getGranularity()
  {    
    return granularity;
  }

  // FIXME?: destructors are NOT called, this is the responsability of the user !
  void freeChunks()
  {
    for (StorageIterator it=storage.begin();it!=storage.end();it++) 
      free(*it);   
    storage.clear();
    allocatedSize.clear();
    freeData=NULL;
    allocEnd=NULL;
    nUsed=0;
    nAllocated=0;   
  }

  void reserve(long n)
  {
    long nLeft=getAllocatedCount()-getUsedCount();    
    if (nLeft<n+1) 
      nReserve = (n+1)-nLeft;
    else 
      nReserve=0;
  }

  void setElementName(std::string elName)
  {
    elementNameStr=elName;
  }

  long getNAllocMin()
  {
    return nAllocMin;
  }

  void setNAllocMin(long nMin)
  {
    nAllocMin=nMin;
  }

  double getAllocFactor()
  {
    return allocFactor;
  }

  void setAllocFactor(double factor)
  {
    allocFactor=factor;
  }
  /*
  void setAllocChunkSize(long allocCount)
  {    
    if (allocCount>1000)
      allocChunkSize = allocCount;
    else
      allocChunkSize = 1000;
      
  } 
  
  void setAllocChunkSizeBytes(long allocSize)
  {
    allocChunkSize = (long)(allocSize/sizeof(T))+1;
  } 
  */
  
  long getUsedCount() const
  {
    return nUsed;
  }

  long getAllocatedCount() const
  {
    return nAllocated;
  }

  long getRecycledCount() const
  {
    return nAllocated-(nUsed+getNSpare());
  }

  long capacity() const
  {
    return nAllocated;
  }

  long size() const
  {
    return nUsed;
  }

  void pop(Type** ptr)
  {
    if (freeData==NULL) newChunk();
    
    T* result = (T*) freeData;
    if (freeData->next == freeData) freeData=NULL;
    else freeData=freeData->next;
  
    ++nUsed;
    new (result) T(); // call the constructor via placement new
    *ptr=result;
    if (result==allocEnd) ++allocEnd;
    //nRecycled--;
  }

protected:
  void recycle_noDestructor(Type* ptr) 
  {    
    MemoryPoolList *tmp = (MemoryPoolList *)ptr;
    tmp->next=freeData;   
    freeData=tmp;
    --nUsed;
  }

public:
  void recycle(Type* ptr) 
  {
    ptr->T::~T();
    recycle_noDestructor(ptr);
    /*
    if (callDestructor) ptr->T::~T(); // call the destructor
    MemoryPoolList *tmp = (MemoryPoolList *)ptr;
    tmp->next=freeData;   
    freeData=tmp;
    --nUsed;

    // if (nRecycled<0) nRecycled=1;
    // else nRecycled++; 
    */ 
  }

  void recycle(Type** ptr) 
  {
    recycle(*ptr);
    (*ptr) = NULL;

    // if (nRecycled<0) nRecycled=1;
    // else nRecycled++;
  }

  void swap(MyType &other)
  {
    std::swap(storage,other.storage);
    std::swap(allocatedSize,other.allocatedSize);
    std::swap(sortedStorage,other.sortedStorage);
    std::swap(sortedCumAllocatedSize,other.sortedCumAllocatedSize);
    std::swap(nUsed,other.nUsed);
    std::swap(nAllocated,other.nAllocated);
    std::swap(allocEnd,other.allocEnd);
    std::swap(nReserve,other.nReserve);
    std::swap(nAllocMin,other.nAllocMin);
    std::swap(allocFactor,other.allocFactor);
    std::swap(freeData,other.freeData);
  }

  // Returns a unique index for each pointer belonging to the pool (recycled or in use)
  // Return -1 if the pointer does not belong to the pool
  // Complexity is ~log(N) with N the number of allocated chunks
  long indexOf(T* ptr) const
  {
    // retrieve the chunk
    typename Storage::const_iterator it = 
      hlp::findFirstHigher(sortedStorage.cbegin(),sortedStorage.cend(),ptr);
    
    if (it==sortedStorage.cbegin()) return -1;
    else --it;

    long storagePos = std::distance(sortedStorage.cbegin(),it);
    long startIndex = (storagePos==0)?0:sortedCumAllocatedSize[storagePos-1];
    long index = startIndex + std::distance(*it,ptr);

    return (index<sortedCumAllocatedSize[storagePos])?index:-1;
  }

  template <class W>
  size_t writeElements(T* el, long N, W* writer, hlp::IsTrue)
  {
    return writer->write(el,N);
  }

  template <class W>
  size_t writeElements(T* el, long N, W* writer, hlp::IsFalse)
  {
    for (long i=0;i<N;++i) T::selfSerialize(&el[i],writer);
    return (N<0)?0:N;
  }

  template <class R>
  size_t readElements(T* el, long N, R* reader, hlp::IsTrue)
  {
    return reader->read(el,N);
  }

  template <class R>
  size_t readElements(T* el, long N, R* reader, hlp::IsFalse)
  {    
    for (long i=0;i<N;++i) T::selfUnSerialize(&el[i],reader);
    return (N<0)?0:N;
  }

  template <class W>
  void serialize(W *writer)
  {    
    writer->writeHeader(classHeader(),classVersion());

    int elSize=sizeof(T);   
    writer->write(&elSize);
  
    writer->write(&elementNameStr);
    writer->write(&nReserve);
    writer->write(&nUsed);
    writer->write(&nAllocated);
    writer->write(&nAllocMin);
    writer->write(&allocFactor);    

    writer->write(&sortedStorage);
    writer->write(&sortedCumAllocatedSize);

    if (nUsed==0) return;

    // nSpare is the number of objects that were never poped
    long nUnused = nAllocated-nUsed;
    long nSpare = std::distance(allocEnd,storage.back()+allocatedSize.back());
    if (nSpare<0) nSpare=0;
    
    writer->write(&nUnused);
    writer->write(&nSpare);
    
    MemoryPoolList *lst = freeData;
    T *ptr = (T*)lst;
    
    // The never poped objects are at the end of the last chunk, so we do not need to
    // write them down
    for (long i=0;i<nUnused;++i)
      {
	if ((ptr<allocEnd)||(ptr>=allocEnd+nSpare)) 
	  writer->write(&ptr);
	lst=lst->next;
	ptr=(T*)lst;
      }
    
    unsigned long prevCumSize=0;    
    typename std::vector<unsigned long>::iterator sz_it=sortedCumAllocatedSize.begin();
    for (StorageIterator it=sortedStorage.begin();it!=sortedStorage.end();++it,++sz_it)
      {	
	unsigned long size = (*sz_it) - prevCumSize;
	prevCumSize = (*sz_it);
	ptr=*it;
	
	// skip never poped objects
	if (ptr==storage.back()) size -= nSpare;	
	
	writer->write(&size);
	writeElements(ptr,size,writer,typename hlp::IsTrueT<IS_POD>::Result());
        //writer->write(ptr,size);
      }   
  }

  template <class R>
  UnserializedPointerUpdater unSerialize(R *reader)
  {
    if (reader == NULL)
      return UnserializedPointerUpdater();
    
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
    
    // Self destruction !
    freeChunks();   
    
    reader->read(&elementNameStr);
    reader->read(&nReserve);
    reader->read(&nUsed);    
    reader->read(&nAllocated);
    reader->read(&nAllocMin);
    double allocFactorFromFile;
    reader->read(&allocFactorFromFile); // use alloc factor from file or the current one ?
        
    reader->read(&sortedStorage);
    reader->read(&sortedCumAllocatedSize);

    // Build the pointer updater so that user can update pointers to pool elements
    UnserializedPointerUpdate *pu = new UnserializedPointerUpdate();

    if (nUsed==0) 
      {
	sortedStorage.clear();
	sortedCumAllocatedSize.clear();
	return UnserializedPointerUpdater(pu);
      }

    // we will reload everything into a single chunk, discarding free unused data
    //T* curChunk = (T*) malloc(sizeof(T)*nUsed); 
    //T* curChunk = (T*)aligned_alloc(ALIGNMENT,sizeof(T)*nUsed);
    T* curChunk = (T*)DICE_ALIGNED_MALLOC(ALIGNMENT,sizeof(T)*nUsed);

    storage.push_back(curChunk);
    allocatedSize.push_back(nUsed);
    allocEnd=curChunk+nUsed; // all elements in the container are used   

    // Use directly 'freeData' in 'pu' so that we do not need to make a copy there later
    std::vector<T*> &puFreeData = pu->getFreeData();
    long nUnused;
    long nSpare;
    reader->read(&nUnused);
    reader->read(&nSpare);
    puFreeData.resize(nUnused-nSpare);
    reader->read(&puFreeData.front(),nUnused-nSpare);     
    std::sort(puFreeData.begin(),puFreeData.end());        
    freeData=NULL;

    // glb::console->print<LOG_STD_ALL>("%s: nUsed=%ld(%ld free), nchunks=%ld\n",
    //  				     elementNameStr.c_str(),nUsed,nUnused,sortedStorage.size());
    
    // Read the chunk into the single allocated new chunk
    // NB: We only restore non freed data so that we get a clean pool !
    // size_t delta=0;     
    // typename std::vector<T*>::iterator free_it = puFreeData.begin();
    // for (StorageIterator it=sortedStorage.begin();it!=sortedStorage.end();++it)
    size_t delta=0;
    typename std::vector<T*>::iterator free_it = puFreeData.begin();
    for (unsigned long index=0;index<sortedStorage.size();++index)
      {	
	unsigned long size;
	reader->read(&size);

	if (index==0) sortedCumAllocatedSize[0]=size;
	else sortedCumAllocatedSize[index]=sortedCumAllocatedSize[index-1]+size;

	// glb::console->print<LOG_STD_ALL>("SZ : cur=%ld, +%ld / %ld\n",
	//  				 delta,size,nUsed);
	// prevCumSize = (*sz_it);
	T* start = sortedStorage[index];
	T* stop = start + size;
	T* cur = start;
	T dump;
	size_t nRead=0;
	
	while (nRead<size)
	  {
	    if (free_it == puFreeData.end()) // the remaining chunks are pristine
	      {
		long dist = std::distance(cur,stop);
		//nRead+=reader->read(curChunk+delta,dist);
		nRead += readElements(curChunk+delta,dist,
				      reader,typename hlp::IsTrueT<IS_POD>::Result());
		delta+=dist;
	      }
	    else // Remaining chunks may be trashed !
	      {
		T* nextFree=*free_it;
		
		if (nextFree>=stop) //read until the end of the chunk
		  {
		    long dist = std::distance(cur,stop);
		    //nRead+=reader->read(curChunk+delta,dist);
		    nRead+=readElements(curChunk+delta,dist,
					reader,typename hlp::IsTrueT<IS_POD>::Result());
		    delta+=dist;		    
		  }
		else if (nextFree==cur) // read one junk data and dump it
		  {
		    //nRead+=reader->read(&dump);
		    nRead+=readElements(&dump,1,
					reader,typename hlp::IsTrueT<IS_POD>::Result());
		    ++free_it;	
		    ++cur;
		    //glb::console->print<LOG_WARNING>("DUMPED 1\n");
		  }
		else // read until next junk data, then read it and dump it
		  {
		    long dist = std::distance(cur,nextFree);
		    //nRead+=reader->read(curChunk+delta,dist);
		    nRead+=readElements(curChunk+delta,dist,
					reader,typename hlp::IsTrueT<IS_POD>::Result());
		    delta+=dist;
		    //nRead+=reader->read(&dump);
		    nRead+=readElements(&dump,1,
					reader,typename hlp::IsTrueT<IS_POD>::Result());
		    ++free_it;
		    cur+=dist+1;
		    //glb::console->print<LOG_WARNING>("DUMPED 1\n");
		  }
	      }	    
	  }

	// if (nRead!=size) 
	//   glb::console->print<LOG_WARNING>("nRead = %ld <= %ld / %ld\n",nRead,size,puFreeData.size());

      }
    nAllocated=nUsed;

    // setup pu correctly (freeData is already OK)
    pu->sortedStorage = sortedStorage;
    pu->sortedCumAllocatedSize = sortedCumAllocatedSize;
    pu->newChunk = curChunk;
    pu->getReady();

    // Check if we failed swapping data bytes if endianness differ
    if (IS_POD)
      {
	bool needByteSwap = reader->getNeedSwap();
	// This happens when swapping is needed and we are using a non primitive type
	if ((!R::template canAutoSwap<T>())&&(needByteSwap))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Structure swapping needed but not implemented yet.\n");
	    exit(-1);
	  }
      }

    // update the sortedStorage and allocatedSize
    sortedStorage=storage;
    sortedCumAllocatedSize=allocatedSize;
    
    return UnserializedPointerUpdater(pu);
  }  

  // FIXME: this is untested !
  // FIXME: write to memory directly ? (use fmemopen ?)
  UnserializedPointerUpdater defrag()
  {
    FILE *tmp = tmpfile(); 

    myIO::BinaryWriterT<> *w = new myIO::BinaryWriterT<>(tmp);
    serialize(w);
    delete w;

    rewind(tmp);

    myIO::BinaryReaderT<> *r = new myIO::BinaryReaderT<>(tmp);
    UnserializedPointerUpdater result = unSerialize(r);
    delete r;

    fclose(tmp);
    return result;
  }
  /*
  long getElementIndex(T *element) const
  {
    // retrieve the chunk
    typename Storage::const_iterator it = 
      hlp::findFirstHigher(sortedStorage.cbegin(),sortedStorage.cend(),element);
    
    if (it==sortedStorage.cbegin()) 
      return -1;
    else 
      --it;

    // get its first and last index    
    long storagePos = std::distance(sortedStorage.cbegin(),it);
    long startIndex = (storagePos==0)?0:sortedCumAllocatedSize[storagePos-1];
    long index = startIndex + std::distance(*it,element);

    return index;
  }
  */
  const std::string& getElementNameStr() {return elementNameStr;}

protected:
  std::string elementNameStr;

  typedef std::vector< T* > Storage; 
  typedef typename Storage::iterator StorageIterator;
  typedef struct MemoryPoolListStruct MemoryPoolList;  
 
  Storage storage; // list of all the chunks
  Storage sortedStorage; // same, but sorted for indexOf

  std::vector<unsigned long> allocatedSize; // size of each chunk  
  std::vector<unsigned long> sortedCumAllocatedSize; // cumulative size of sorted chunks 

  T* allocEnd; // A pointer to the first never poped adress
  long nReserve; // how many at least should be allocated next time newChunk is called
  long nUsed; // how many have been poped
  long nAllocated; // how many have been allocated in total
  //long nRecycled; // how many pointers have been recycled (can be negative)
  long nAllocMin; // how many at least should be allocated whenever newChunk is called
  double allocFactor; // what fraction of the current size should be allocated when newChunk is called
  MemoryPoolList *freeData; // list of free objects
  long granularity; // The number of object allocated per page must be a multiple of granularity

  long getNSpare() const
  {
    long nSpare;
    if (storage.size() == 0)
      nSpare=0;
    else
      nSpare = std::distance(allocEnd,
			     storage.back()+allocatedSize.back());
    if (nSpare<0) nSpare=0;

    return nSpare;
  }
 
  virtual void newChunk(size_t forcedSize=0)
  {  
    long allocChunkSize = getAllocatedCount()*allocFactor +1;

    long nAlloc=(nReserve>allocChunkSize) ? nReserve : allocChunkSize;
    if (nAlloc<nAllocMin) nAlloc = nAllocMin;

    long complement=granularity-(nAlloc%granularity);
    nAlloc += complement;

    if (forcedSize>0) nAlloc = forcedSize;

    //printf("Alloc : %ld (%ld,%g)\n",nAlloc,getAllocatedCount(),allocFactor);
    //glb::console->printNewLine<LOG_DEBUG>("%s: (%d) ? %ld : %ld\n",elementNameStr.c_str(),(int)(nReserve>allocChunkSize),(long)nReserve,(long)allocChunkSize);

    // use malloc NOT new, as we don't want to call constructors at all
    // or destructors when using delete ...
    //T* curChunk = (T*) malloc(sizeof(T)*nAlloc); 
    //T* curChunk = (T*)aligned_alloc(ALIGNMENT, sizeof(T)*nAlloc);
    T* curChunk = (T*)DICE_ALIGNED_MALLOC(ALIGNMENT, sizeof(T)*nAlloc);

    if (curChunk == NULL)
      {
	if (forcedSize==0) 
	  {
	    glb::console->print<LOG_WARNING>("Pool '%s': Could not allocate %ld bytes.",
					     elementNameStr.c_str(),
					     sizeof(T)*nAlloc);
	    glb::console->print<LOG_WARNING>("  Will try allocating %ld bytes ...\n",
					     (sizeof(T)*nAlloc)/2);
	    newChunk(nAlloc/2);
	    return;
	  }
	
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Pool '%s': Could not allocate %ld bytes.",
				       elementNameStr.c_str(),
				       sizeof(T)*nAlloc);
	exit(-1);
      }

    MemoryPoolList *cur=(MemoryPoolList*) &(curChunk[0]);
    for (long i=1;i<nAlloc;i++)
      {
	MemoryPoolList *next = (MemoryPoolList*) &(curChunk[i]);
	cur->next=next;
	cur=next;
      }
    cur->next=cur; // marks end of list

    if (freeData==NULL) freeData = (MemoryPoolList*) &(curChunk[0]);
    else
      {
	cur->next=freeData;
	freeData=(MemoryPoolList*) &(curChunk[0]);
      }

    storage.push_back(curChunk);
    allocatedSize.push_back(nAlloc);
    
    // update the sorted chunck position and sizes
    if (sortedStorage.size()>0)
      {
	typename Storage::iterator it = 
	  hlp::findFirstHigher(sortedStorage.begin(),sortedStorage.end(),storage.back());
	
	long index = std::distance(sortedStorage.begin(),it);
	
	// this is used for indexOf and serialization
	sortedStorage.insert(it,storage.back());
	sortedCumAllocatedSize.insert(sortedCumAllocatedSize.begin()+index,
				      allocatedSize.back());
    
	if (index!=0) 
	  sortedCumAllocatedSize[index]+=sortedCumAllocatedSize[index-1];

	for (long i=index+1;i<sortedCumAllocatedSize.size();++i)
	  sortedCumAllocatedSize[i]+=allocatedSize.back();
      }
    else
      {
	sortedStorage=storage;
	sortedCumAllocatedSize=allocatedSize;
      }
    //glb::console->flushBuffer<LOG_STD_ALL>();
    // for (long i=0;i<sortedStorage.size();++i)
    //   {
    // 	glb::console->print<LOG_STD_ALL>(" (%s:%ld - %ldB)",elementNameStr.c_str(),
    // 						 (long)sortedStorage[i],
    // 						 (long)sortedCumAllocatedSize[i]);
    //   }
    //glb::console->printToBuffer<LOG_STD_ALL>("\n");
    //glb::console->flushBuffer<LOG_STD_ALL>();

    nAllocated+=nAlloc;
    allocEnd=curChunk;
    nReserve=0;
    
    glb::console->printNewLine<LOG_DEBUG>("MemoryPool: allocated %ld %s of size %ld. (%g Go)\n",
				     nAlloc,elementNameStr.c_str(),sizeof(T),
				     (double)(nAlloc*sizeof(T))/(1l<<30));
    //printf("Allocated %ld el from %ld to %ld\n",curChunkPtr->size(),(long)&(curChunk.front()),(long)(allocEnd+curChunkPtr->size()));
  }
};


// Definition of UnserializedPointerUpdate
// Used to update pointers to memory pool elements after unserializing
template <class T, bool IS_POD>
class MemoryPoolT<T,IS_POD>::UnserializedPointerUpdate
{
public:
  typedef T Data;
  typedef T* DataPtr;
  friend class MemoryPoolT<T,IS_POD>;
            
  template <class RT>
  RT* operator()(RT* ptr) const
  {
    return static_cast<RT*>(convert(static_cast<DataPtr>(ptr)));
  }
    
private:
  /*
  PointerUpdate(const std::vector<T*> &sortedStorage_, 
		const std::vector<unsigned long> &sortedCumAllocatedSize_,
		const DataPtr newChunk_):
    sortedStorage(sortedStorage_),
    sortedCumAllocatedSize(sortedCumAllocatedSize_),
    newChunk(newChunk_)
  {}
  */
  UnserializedPointerUpdate()
  {
    sortedStorage.clear();
    sortedCumAllocatedSize.clear();
    newChunk=NULL;
  }
  /*
  void set(const std::vector<T*> &sortedStorage_, 
	   const std::vector<unsigned long> &sortedCumAllocatedSize_,
	   const std::vector<T*> freeData_,
	   const DataPtr newChunk_)
  {
    sortedStorage=sortedStorage_;
    sortedCumAllocatedSize=sortedCumAllocatedSize_;
    freeData=freeData_;
    newChunk=newChunk_;
    getReady();
  }
  */
  std::vector<T*> &getFreeData() {return freeData;}
  
  void getReady()
  {
    if (freeData.size()==0) 
      freeData.push_back(sortedStorage.back()+sortedCumAllocatedSize.back()+1);
  }
  
  DataPtr convert(const DataPtr ptr) const
  {
    // retrieve the chunk
    typename Storage::const_iterator it = 
      hlp::findFirstHigher(sortedStorage.cbegin(),sortedStorage.cend(),ptr);
    
    if (it==sortedStorage.cbegin()) 
      return static_cast<DataPtr>(NULL);
    else 
      --it;

    // get its first and last index    
    long storagePos = std::distance(sortedStorage.cbegin(),it);
    size_t startIndex = (storagePos==0)?0:sortedCumAllocatedSize[storagePos-1];
    size_t index = startIndex + std::distance(*it,ptr);

    // glb::console->print<LOG_STD_ALL>
    //   ("Updating pointer @%ld: found @chunk %ld, index= %ld < %ld < %ld \n",
    //    (long)ptr,storagePos,
    //    startIndex,index,sortedCumAllocatedSize[storagePos]);
       
    if (index<sortedCumAllocatedSize[storagePos])
      {
	if (ptr<freeData.front()) 
	  return newChunk+index;
	else if (ptr>freeData.back())
	  return newChunk+index-freeData.size();
	else
	  {
	    typename Storage::const_iterator it = 
	      hlp::findFirstHigher(freeData.cbegin(),freeData.cend(),ptr);
	    index -= std::distance(freeData.cbegin(),it);
	    --it;
	    
	    if (ptr!=(*it))
	      return newChunk+index;
	    else 
	      return static_cast<DataPtr>(NULL);
	  }
      }
    else return static_cast<DataPtr>(NULL);
  }
  DataPtr newChunk;
  
  // Each vector must be sorted in ASCENDING order !
  std::vector<T*> sortedStorage;
  std::vector<unsigned long> sortedCumAllocatedSize;
  std::vector<T*> freeData;
};
/** \}*/
#include "../../internal/namespace.footer"
#endif
