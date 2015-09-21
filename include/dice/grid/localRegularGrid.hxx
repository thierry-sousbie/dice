#ifndef __LOCAL_REGULAR_GRID_HXX__
#define __LOCAL_REGULAR_GRID_HXX__

#include <map>
#include <vector>
#include <stdio.h>

#include "./valLocationType.hxx"
#include "./scale.hxx"
#include "./regularGridNavigation.hxx"
#include "./localRegularGridParams.hxx"
#include "./regularGridFieldLayout.hxx"
#include "./gridKernels.hxx"
#include "./regularGridSymmetryEnum.hxx"

#include "../geometry/boundaryType.hxx"
#include "../geometry/geometricProperties.hxx"

#include "../IO/vtkRectilinearGridWriter.hxx"

#include "../finiteElements/gaussQuadrature.hxx"

#include "./internal/subgridIterator.hxx"
//#include "./internal/interpolateCIC.hxx"
#include "./internal/lrg_dispatchers.hxx"
#include "./internal/lrg_quadratureFunctorAdapter.hxx"
#include "./internal/symmetryFunctors.hxx"

/**
 * @file 
 * @brief A regular grid class
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

/**
 * \class LocalRegularGridT
 * \brief  A regular grid class
 * \tparam ND The number of dimensions
 * \tparam DT The type of data to store at grid sites
 * \tparam BT The boudary conditions type ( see BoundaryType )
 * \tparam FML defines the memory layout of different fields in the grid ( see 
 * regularGridFieldLayout )
 * \note Data is stored in column-major order (i.e. Fortran order) : the first dimension
 * index changes first when traversing the array linearly in memory.
 */
template <long ND, 
	  typename DT=double,
	  long BT = BoundaryType::NONE,
	  long FML = regularGridFieldLayout::CONSECUTIVE>
class LocalRegularGridT {
public:
  typedef LocalRegularGridT<ND,DT,BT,FML> MyType;  

  static const long BOUNDARY_TYPE = BT;
  static const int  NDIM = ND;

  static const int  IS_INTERLEAVED = (FML==regularGridFieldLayout::INTERLEAVED);
  static const int  IS_PERIODIC = (BT==BoundaryType::PERIODIC);
  static const int  LAYOUT = FML;
 
  static std::string classHeader() {return "local_regular_grid";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  typedef DT value_type;
  typedef DT Data;

  typedef RegularGridNavigation GridNav;
  typedef typename RegularGridNavigation::Direction Direction;

  friend class internal::SubgridIteratorT<MyType>;
  template <int order, int diff_order> friend struct internal::KernelDispacherT;

  // CONSECUTIVE : iterate over each stored elements, as stored in memory
  // INTERLEAVED : iterate of each elements nfield by nfields, as stored in memory
  typedef internal::SubgridIteratorT<MyType> iterator;

  // iterate over each field of each element in order (E0_0,..,E0_F,E1_0,..,E1_F,....,EN_F)
  // This is optimized for INTERLEAVED, slow for CONSECUTIVE
  typedef internal::field_iterator<MyType> fieldIterator;

  // iterate over the values of one field in order (E0_k,E1_k,E2_k,....,EN_k)  
  typedef internal::oneField_iterator<MyType> oneFieldIterator;

  // iterate the same way data is stored in memory 
  typedef internal::inOrder_iterator<MyType> inOrderIterator;
 
  typedef ScaleT<double> Scale; 
  
  typedef typename Scale::ScaleType  ScaleType;
  typedef typename Scale::ScaleTypeV ScaleTypeV;    
  //typedef ValLocationType ValLocationType;
  //typedef ValLocationTypeV ValLocationTypeV;

  typedef LocalRegularGridParamsT<ND> Params;
  
protected:
  Params params; 
  long nVertices;
  long nCells;
  long nValues;
  long nFields;
  
  long nAllocated;
  
  double xMin[NDIM];
  double xMax[NDIM];
  double boxSize[NDIM];
  double boxSizeInv[NDIM];
  double arrDim[NDIM];

  std::vector<double> vertexCoord[NDIM];
  std::vector<double> cellCoord[NDIM];
  std::vector<double> valueCoord[NDIM];
  std::vector<double> vertexDelta[NDIM];
  std::vector<double> cellDelta[NDIM];
  std::vector<double> valueDelta[NDIM];

  double constantVolumeElement;

  long vertexStride[NDIM+1];
  long cellStride[NDIM+1];
  long valueStride[NDIM+1];

  // Kernels
  mutable gridKernel::NGP<ND> gridKernel_NGP;
  mutable gridKernel::CIC<ND> gridKernel_CIC;
  mutable gridKernel::TSC<ND> gridKernel_TSC;
  
  mutable gridKernel::InterpolateCentralDiffT<ND,1> gridKernel_grad_CIC;
  mutable gridKernel::InterpolateCentralDiffT<ND,2> gridKernel_grad_TSC;  

  std::vector<long> integrationPoints[1<<NDIM];
  
  std::string dataName;

  Data *arr;
  bool ownArr;
  
private:
  bool initialized;
  typedef std::map<std::string,double> valueMapT;
  typedef typename valueMapT::iterator valueMapItT;
  valueMapT valueMap;

public:

  bool registerValue(const char *name, double value, bool replace=true)
  {
    std::string pname(name);
    std::pair<valueMapItT,bool> res=
      valueMap.insert(std::make_pair(pname,value));

    if ((!res.second)&&(replace)) {
      valueMap.erase(res.first);
      res=valueMap.insert(std::make_pair(pname,value));
      return false;
    }
    
    return true;
  }

  bool getRegisteredValue(const char *name, double &value)
  {
    std::string pname(name);
    valueMapItT res=valueMap.find(pname);
    if (res==valueMap.end()) return false;
    value=res->second;
    return true;
  }

  /** \brief Save the grid to a VTK RectilinearGrid file
   *  \param fname The name of the file, or NULL for default name
   */
  void toVtk(const char *fname=NULL)
  {
    IO::VtkRectilinearGridWriterT<MyType> vtkWriter(this,fname);
    vtkWriter.write();
  }

  std::string getName() const
  {
    return dataName;
  }

  void setName(const char *name)
  {
    dataName=name;
  }

  template <class K>
  void initializeInterpolationKernel(K &kernel, bool setSpacing=true) const
  {
    kernel.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
    
    if (setSpacing)
      {
	for (int i=0;i<NDIM;++i)
	  kernel.setSpacing(getValueDelta(i)[0],i);
      }
  }

private:  

  void initializePrivateKernels()
  {
    initializeInterpolationKernel(gridKernel_NGP);
    initializeInterpolationKernel(gridKernel_CIC);
    initializeInterpolationKernel(gridKernel_TSC);
    initializeInterpolationKernel(gridKernel_grad_CIC);
    initializeInterpolationKernel(gridKernel_grad_TSC);
    // gridKernel_NGP.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
    // gridKernel_CIC.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
    // gridKernel_TSC.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
    // gridKernel_grad_CIC.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
    // gridKernel_grad_TSC.initialize(params.resolution,valueStride,(IS_INTERLEAVED)?getNFields():1);
  }

  void initialize(value_type *myArr=NULL, long minElementsCount=0)
  {    
    int i,j;
    /*
    if (initialized) 
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("Trying to initialize LocalRegularGrid twice !\n");
      }
    */
    nVertices=1; 
    nCells=1; 
    nValues=1;
    nFields=params.nFields;
        
    vertexStride[0]=1;//1*nFields;
    cellStride[0]=1;//1*nFields;
    valueStride[0]=1;//1*nFields;

    constantVolumeElement = 1.0;
  
    for (i=0;i<NDIM;i++) {
      xMin[i]=params.x0[i];
      xMax[i]=params.x0[i]+params.delta[i];
      vertexCoord[i]=Scale::genScale
	(xMin[i],xMax[i],params.resolution[i],params.scale[i],ValLocationTypeV::VERTEX);
      cellCoord[i]=Scale::genScale
	(xMin[i],xMax[i],params.resolution[i],params.scale[i],ValLocationTypeV::CELL);
      valueCoord[i]=Scale::genScale
	(xMin[i],xMax[i],params.resolution[i],params.scale[i],params.valLocation[i]);

      if ((params.valLocation[i]==ValLocationTypeV::VERTEX)&&(IS_PERIODIC))
	  valueCoord[i].resize(valueCoord[i].size()-1);

      vertexDelta[i]=Scale::genScaleDelta
	(xMin[i],xMax[i],params.resolution[i],params.scale[i],ValLocationTypeV::VERTEX);
      cellDelta[i]=Scale::genScaleDelta
	(xMin[i],xMax[i],params.resolution[i],params.scale[i],ValLocationTypeV::CELL);
      cellDelta[i].push_back(0);

      if (params.valLocation[i]==ValLocationTypeV::CELL)
	valueDelta[i] = cellDelta[i];
      else
	valueDelta[i] = vertexDelta[i];

      constantVolumeElement*=valueDelta[i][0];

      nVertices*=vertexCoord[i].size();
      nCells*=cellCoord[i].size();
      nValues*=valueCoord[i].size();

      cellStride[i+1]=cellStride[i]*cellCoord[i].size();
      vertexStride[i+1]=vertexStride[i]*vertexCoord[i].size();
      valueStride[i+1]=valueStride[i]*valueCoord[i].size();      
    }   

    nAllocated = std::max(std::max(minElementsCount,params.minElementsCount),
			  nValues*nFields);
   
    if ((arr!=NULL)&&(ownArr)) free(arr);
    
    if (myArr!=NULL)
      {
	ownArr=false;
	arr=myArr;
      }
    else
      {
	ownArr=true;
	arr=(value_type*)calloc(nAllocated,sizeof(value_type));
	int res=(arr==NULL);
	//int res=posix_memalign((void**)&arr,64,nAllocated*sizeof(value_type));
	if (res!=0)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>
	      ("Trying to allocate %ld bytes of 64 bytes aligned memory.\n",nAllocated*sizeof(value_type));
	    glb::console->print<LOG_ERROR>
	      ("posix_memalign returned error %d.\n",res);
	    exit(-1);
	  }

	// FIXME: For very large arrays (e.g. 64Gb), this seems to take forever on some systems with posix_memalign. i don't know why, have to investigate ...
	//memset(arr,0,nAllocated*sizeof(value_type));
      }

    for (int k=0;k<(1<<NDIM);k++) integrationPoints[k].clear();
    for (int k=0;k<(1<<NDIM);k++)
      {	
	for (i=0;i<(1<<NDIM);i++)
	  {	
	    long ptr=0;	    
	    for (j=0;j<NDIM;j++)
	      {
		if (i&(1<<j))
		  {
		    if (k&(1<<j)) break;
		    ptr+=valueStride[j];
		  }
	      }
	    if (j==NDIM) integrationPoints[k].push_back(ptr*nFields);
	  }
      }
    
    for (i=0;i<NDIM;i++)
      {
	arrDim[i]=getValueCoord(i).size();
	boxSize[i]=getVertexCoord(i).back()-getVertexCoord(i).front();	
	boxSizeInv[i]=1.0L/boxSize[i];
      }

    initializePrivateKernels();

    initialized=true;
  }

public:
  value_type *getDataPtr() {return arr;}

  value_type *getDataPtr(long i)
  {
    if (IS_INTERLEAVED)
      return &arr[i*nFields];
    else
      return &arr[i];
  }
  
  value_type *getDataPtr(long i, int field) 
  {
    if (IS_INTERLEAVED)
      return &arr[i*nFields+field];
    else
      return &arr[i+field*getNValues()];
  }

  template <typename T>
  value_type *getDataPtr(T* w) 
  {
    long index=w[0]*valueStride[0];
    for (int i=1;i<NDIM;i++) 
      index+=valueStride[i]*w[i];
     
    return getDataPtr(index);   
  }

  template <typename T>
  value_type *getDataPtr(T* w, int field) 
  {
    long index=w[0]*valueStride[0];

    for (int i=1;i<NDIM;i++) 
      index+=valueStride[i]*w[i];
   
    return getDataPtr(index,field);   
  }

  template <typename TD, typename TI>
  void coordsToICoords(TD *coords, TI *iCoords)
  {
    for (int i=0;i<NDIM;++i)
      iCoords[i] = ((coords[i]-getOrigin(i)) * getBoxSizeInv(i) * getValueCoord(i).size());
  }

  template <typename TI,typename TD>
  void iCoordsToCoords(TI *iCoords, TD *coords)
  {
    for (int i=0;i<NDIM;++i)
      {
	coords[i]=getValueCoord(i)[iCoords[i]];      
      }
  }

  long getCellStride(int i) const {return cellStride[i];}
  long getVertexStride(int i) const {return vertexStride[i];}
  long getValueStride(int i) const {return valueStride[i];}
  const long *getValueStride() const {return valueStride;}
  const std::vector<long> &getIntegrationPoints(int id=0) const 
  {return integrationPoints[id];}

  virtual long getLowMargin(int i) const {return params.lowMargin[i];}
  virtual long getHighMargin(int i) const {return params.highMargin[i];}

  long getNVertices() const {return nVertices;}
  long getNCells() const {return nCells;}
  long getNValues() const {return nValues;}
  long getNFields() const {return nFields;}
  long getArrayDim(int i) const {return arrDim[i];}
  double getBoxSize(int i) const {return boxSize[i];}
  double getBoxSizeInv(int i) const {return boxSizeInv[i];}
  const double *getBoxSizeInv() const {return boxSizeInv;}

  const std::vector<double>& getCellCoord(long i) const 
  {return cellCoord[i];}
  const std::vector<double>& getVertexCoord(long i) const 
  {return vertexCoord[i];}
  const std::vector<double>& getValueCoord(long i) const 
  {return valueCoord[i];}  

  const std::vector<double>& getCellDelta(long i) const 
  {return cellDelta[i];}
  const std::vector<double>& getVertexDelta(long i) const 
  {return vertexDelta[i];}
  const std::vector<double>& getValueDelta(long i) const 
  {return valueDelta[i];}

  double getConstantVolumeElement() const
  {
    return constantVolumeElement;
  }

  ValLocationType getValLocation(int dim) const {return params.valLocation[dim];}

  const Params &getParams() const {return params;} 
  bool isInitialized() const {return initialized;}
  
  int getPosition(int dim) const
  {
    if (!params.haveParentGrid)
      return 0;
    else
      return params.position[dim];
  }

  int getResolution(int dim) const
  {
    return params.resolution[dim];
  }

  double getSize(int dim) const 
  {
    return params.delta[dim];
  }

  double getOrigin(int dim) const 
  {
    return params.x0[dim];
  }

  const double *getOrigin() const 
  {
    return params.x0;
  }

  /*
  double getSize(int dim) const
  {
    return valueCoord[dim].size();
  }
  */
  int getGlobalResolution(int dim) const
  {
    if (!params.haveParentGrid)
      return getResolution(dim);
    else
      return params.parentResolution[dim];
  }

  template <class T>
  void getBoundingBox(T bbox[NDIM][2]) const
  {
    for (int i=0;i<NDIM;++i)
      {
	bbox[i][0] = params.x0[i];
	bbox[i][1] = params.x0[i] + params.delta[i];
      }
  }

  template <class T>
  void getGlobalBoundingBox(T bbox[NDIM][2]) const
  {
    if (!params.haveParentGrid)
      getBoundingBox(bbox);
    else
      {
	for (int i=0;i<NDIM;++i)
	  {
	    bbox[i][0] = params.parentX0[i];
	    bbox[i][1] = params.parentX0[i] + params.parentDelta[i];
	  }
      }
  }

  void reinterpretNFields(int N=0)
  {
    if (N<=0)
      nFields=params.nFields;
    else nFields=N;
    initializePrivateKernels();
  }
  
public:
  
  explicit LocalRegularGridT(const char *dataName_="value"):
    dataName(dataName_),
    arr(NULL),
    ownArr(false),
    initialized(false)
  {
    
  }

  virtual ~LocalRegularGridT() 
  {
    if ((arr!=NULL)&&(ownArr)) free(arr);
  }

  // This builds a clone and data is always copied to the cloned verion (pointer to
  // data are NOT shared)
  void clone(MyType &cloned, bool cloneExtraElements=true) const
  {
    if ((cloned.arr != NULL)&&(cloned.ownArr))
      free(cloned.arr);  
    
    if (initialized)
      {
	cloned = *this;
	long nClonedElements=nAllocated;
	if (!cloneExtraElements) 
	  {
	    nClonedElements=nValues*nFields;
	    cloned.params.minElementsCount=nClonedElements;	    
	  }
	cloned.arr=NULL; // if we don't do this arr will be freed again !!!!!
	//if (cloned.ownArr)
	  {
	    cloned.initialize(NULL,nClonedElements);
	    std::copy(arr, arr+nClonedElements, cloned.arr);
	  }
      }
    else
      {
	cloned.arr=NULL;
	cloned.initialized=false;
	cloned.ownArr=false;
	cloned.dataName=dataName;
      }
    
  }

  template <class Source>
  long copyRawBuffer(Source &source)
  {
    long count = std::min(nAllocated,source.nAllocated);
    std::copy(source.arr,source.arr+count,arr);
    return count;    
  }

  //private:
  
  LocalRegularGridT( const MyType& other ):
    dataName("copyConstructed"),
    arr(NULL),
    ownArr(false),
    initialized(false)
  {
    other.clone(*this);
  }

private:
  // Prevent public copy but keep default copy for clone implementation
  LocalRegularGridT& operator=(const MyType&) = default;

protected:

  template <class W>
  void writeInit(W *writer)
  {    
    params.write(writer); 
    writer->write(&nAllocated);
    writer->write(&nFields);
    writer->write(dataName);
  }

  template <class W>
  void writeData(W *writer)
  {
    writer->writeHeader(classHeader(),classVersion());
    //writer->write(arr,nValues*nFields);
    writer->write(arr,nAllocated);
  }

public:

  template <class R>
  void initialize(R* reader)
  { 
    Params p;

    p.read(reader);
    reader->read(&nAllocated);
    initialize(p,nAllocated);
    reader->read(&nFields);    
    reader->read(dataName);
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);       
    //reader->read(arr,nValues*nFields);
    reader->read(arr,nAllocated);
  }

  void initialize(const Params &params_, value_type* myArr=NULL)
  {
    params=params_;   
    initialize(myArr);
  }

  template <class W>
  void write(W *writer)
  {
    writeInit(writer);
    writeData(writer);
  }
  
  iterator begin(bool includeMargin=false)
  {
    if (includeMargin)
      return iterator(this,GridNav::undefined());      
    else
      return iterator(this,GridNav::dir());
  }

  iterator end(bool includeMargin=false)
  {
    if (includeMargin)
      return iterator(this,GridNav::undefined(),true);    
    else
      return iterator(this,GridNav::dir(),true);
  }

  iterator margin_begin(Direction dir=GridNav::dir())
  {
    return iterator(this,dir);
  }

  iterator margin_end(Direction dir=GridNav::dir())
  {
    return iterator(this,dir,true);
  }

  long margin_size(Direction dir=GridNav::dir())
  {
    return iterator::boxSize(margin_begin(dir));
  }

  iterator subbox_begin(const int *iMin, const int *iMax)
  {
    return iterator(this,iMin,iMax);
  }

  iterator subbox_end(const int *iMin, const int *iMax)
  {
    return iterator(this,iMin,iMax,true);
  } 

  long subbox_size(const int *iMin,const int *iMax)
  {    
    return iterator::boxSize(subbox_begin(iMin,iMax));
  }

  iterator innerMargin_begin(const int *delta, Direction dir=GridNav::dir())
  {
    int iMin[NDIM];
    int iMax[NDIM];
    int i;
  
    for (i=0;i<NDIM;i++)
      {
	int d=delta[i];
	if (d==0) d=valueCoord[i].size()-params.highMargin[i]-params.lowMargin[i];
	else if (d<0) d=valueCoord[i].size();

	if (dir&GridNav::dir(i,1)) 
	  {
	    iMax[i]=valueCoord[i].size()-params.highMargin[i];
	    iMin[i]=iMax[i] - d;
	  }
	else
	  {
	    iMin[i]=params.lowMargin[i];
	    iMax[i]=iMin[i]+d;
	  }
      }
    
    return subbox_begin(iMin,iMax);
  }

  iterator innerMargin_end(const int *delta, Direction dir=GridNav::dir())
  {
    return iterator(this,dir,true);
  }

  long innerMargin_size(const int *delta, Direction dir=GridNav::dir())
  {
    return iterator::boxSize(innerMargin_begin(delta,dir));
  }
  
  template <class otherItT>
  iterator importIterator(const otherItT &it)
  {   
    if (params.haveParentGrid)
      return iterator(this,it,params.parent_which);
    else
      return iterator(this,it);
  }

  int getParentNDim() {return params.parentNDim;}
  
  template <class itT>
  void interleave(const itT start,const itT stop, double fac=1.)
  {    
    if (std::distance(start,stop)!=nFields)
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("nFields(%d) and number of grids(%ld) do not match.\n",nFields,std::distance(start,stop));
	exit(-1);
      }

    long i=0;
    value_type *other[nFields];

    for (itT it=start;it!=stop;it++)
      other[i++]=&(*it);

    return interleave(other,&other[nFields],fac);
  }

  void interleave(MyType* start, MyType* stop, double fac=1.)
  {
    if (std::distance(start,stop)!=nFields)
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("nFields(%d) and number of grids(%ld) do not match.\n",nFields,std::distance(start,stop));
	exit(-1);
      }    
    
    iterator it=begin();
#pragma omp parallel for
    for (long s=0; s<nFields;s++)
      {
	oneFieldIterator it(begin(),s);
	const oneFieldIterator it_end(end(),s);
	iterator oit=start[s].importIterator(it);
	
	for (;it!=it_end;it++,oit++)
	  (*it)=fac*(*oit);
      }
  }

  // FIXME: openMP
  void getMinMax(double &min, double &max)
  {
    inOrderIterator it = inOrderIterator(begin());
    const inOrderIterator it_end = inOrderIterator(end());
    
    min=(*it);
    max=(*it);

    for (;it!=it_end;++it)
      {
	if (min>(*it)) min=(*it);
	if (max<(*it)) max=(*it);
      }   
  }  

  // FIXME: openMP
  void subMultiply(double sub, double mult)
  {
    inOrderIterator it = inOrderIterator(begin());
    const inOrderIterator it_end = inOrderIterator(end());

    for (;it!=it_end;++it)
      {
	(*it) = ((*it)-sub)*mult;
      }   
  }

  template <class Functor>
  void visit(const Functor &f, int field=0, int nThreads=glb::num_omp_threads)
  {
    
#pragma omp parallel num_threads(nThreads)
      {	
	int th=omp_get_thread_num();
	oneFieldIterator it(begin(),field);
	const oneFieldIterator it_end(end(),field);
	// inOrderIterator it(begin());
	// const inOrderIterator it_end(end());
	it+=th;	
	for (;it!=it_end;it+=nThreads)
	  f(it,th);
      }           
   
  }

  /*
  void multiply(double factor)
  {
#pragma omp parallel for 
    for (long i=0;i<nValues;i++)
      {
	arr[i]*=factor;
      }
  }
  

  void copyRaw(value_type *d, long N)
  {
    if (params.minStorageSize<N)
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("insufficient space to copy, please call set setMinStorage ...\n");
	exit(-1);
      }
    std::copy(d,d+N,arr);
  }

  template <class otherGridT>
  void copyRaw(otherGridT &g)
  {
    copyRaw(g.getDataPtr(),g.getNVal());
  }
  */
  void erase()
  {
    memset(arr,0,nAllocated*sizeof(value_type));
  }

  void reserve(long N)
  {
    long nMin = nValues*nFields;
    long newSize = N;    
    if (newSize<nMin) newSize=nMin;
    
    if (newSize != nAllocated)
      arr=(value_type*)realloc(arr,sizeof(value_type)*newSize);

    nAllocated = newSize;
    /*
    long after=std::max(N,nValues*nFields);
    long before=std::max(parms.minElementsCount,nValues*nFields);
    
    if (before==after) return;
    //printf("Storage before: (%ld,%ld)\n",params.minStorageSize,before);
    if (after<=nValues*nFields)
      params.minElementsCount=nValues*nFields;
    else 
      params.minElementsCount=after;

    if (params.minElementsCount > before)
      arr=(value_type*)realloc(arr,sizeof(value_type)*after);
    */
  }

  template <class AMR, class I, class T>
  void addVoxel(I index, long level, T value,
		long gridLevel,
		long gridPosition[NDIM])
  {
    typedef typename AMR::ICoord ICoord;    
    
    ICoord iCoord[NDIM];
    long min[NDIM];
    long max[NDIM];
    long delta[NDIM];
    Data *ptr=arr;
    bool skip=false;

    long voxelILen=(gridLevel>=level)?(1L<<(gridLevel-level)):1;//-(1L<<(level-gridLevel));
    long voxelHalfILen=(voxelILen>>1);//(voxelILen>=0)?(voxelILen>>1):-((-voxelILen)>>1);
    AMR::index2ICoords(index,iCoord,gridLevel);    
    for (int i=0;i<NDIM;++i)
      {
	long c=long(iCoord[i])-voxelHalfILen-gridPosition[i];
	max[i]=c+voxelILen;
	min[i]=c;	
	if (min[i]<0)
	  {
	    min[i]=0;
	    if (min[i]>=max[i]) skip=true;
	  }
	if (max[i]>params.resolution[i]) 
	  {
	    max[i]=params.resolution[i];
	    if (min[i]>=max[i]) skip=true;
	  }
	delta[i]=max[i]-min[i];
	ptr+=min[i]*valueStride[i];
      }

    
    /*
    glb::console->print<LOG_STD_ALL>("Printing voxel @level%ld for grid @%ld: i=[%ld %ld %ld] P=[%ld %ld %ld] [%ld,%ld][%ld,%ld][%ld,%ld]\n",level,gridLevel,iCoord[0],iCoord[1],iCoord[2],gridPosition[0],gridPosition[1],gridPosition[2],min[0],max[0],min[1],max[1],min[2],max[2]);
    //return;
    */
    //Data value = static_cast<Data>(val);
    /*
    if ((std::distance(arr,ptr)==(113818-113664))&&(glb::mpiComWorld->rank()==3))
      printf("Adding [%ld %ld] [%ld %ld] [%ld %ld] (%g)+=%g -> %g\n",
	     min[0],min[1],max[0],max[1],delta[0],delta[1],
	     *ptr,value,(*ptr)+value);
    */

    if (!skip) 
      {
	if (voxelILen==1) (*ptr)+=value;
	else
	  {	
	    if (NDIM==1)
	      {	    
		for (long i=0;i<delta[0];++i)
		  ptr[i]+=value;
	      }
	    else if (NDIM==2)
	      {
		long dec=0;
		const long d1=valueStride[1];
		for (long j=0;j<delta[1];++j)
		  {		
		    for (long i=dec;i<dec+delta[0];++i)
		      {
			ptr[i]+=value;
		      }
		    dec+=d1;	
		  }
	      }
	    else
	      {
		long dec=0;
		const long d1=valueStride[1];
		const long d2=valueStride[2]-delta[1]*valueStride[1];
		for (long k=0;k<delta[2];++k)
		  {		
		    for (long j=0;j<delta[1];++j)
		      {		    
			for (long i=dec;i<dec+delta[0];++i)
			  {
			    //printf("%ld -> +%g |",i,value);
			    ptr[i]+=value;
			  }
			dec+=d1;
		      }
		    dec+=d2;
		  }
		//printf("\n");
	      }
	  }
      }// skip

  }

  template <class AMR>
  void addRemoteVoxels(const std::vector<typename AMR::Voxel::MpiExchangeData> &voxels,
		       long gridLevel,
		       long gridPosition[NDIM], 
		       int nThreads=glb::num_omp_threads)
  {
    addRemoteVoxels<AMR>(&voxels[0],voxels.size(),gridLevel,gridPosition,nThreads);
  }

  template <class AMR>
  void addRemoteVoxels(const typename AMR::Voxel::MpiExchangeData *voxels, 
		       unsigned long count,
		       long gridLevel,
		       long gridPosition[NDIM], 
		       int nThreads=glb::num_omp_threads)
  {            
    //glb::console->printFlush<LOG_STD>("Adding %ld remote voxels @[%ld %ld]\n",count,gridPosition[0],gridPosition[1]);
    //unsigned long ompDelta=std::max(1UL,count/nThreads/10);

#pragma omp parallel for num_threads(nThreads) //schedule(dynamic,ompDelta)
    for (long i=0;i<count;++i)
      {	
	addVoxel<AMR>(voxels[i].index,voxels[i].level,voxels[i].data,
		      gridLevel,gridPosition);
      }    
    //glb::console->printFlush<LOG_STD>("Adding done.\n");
  }
  
  template <class AMR>
  void addVoxels(typename AMR::Voxel * const *voxels,
		 unsigned long count,
		 long gridLevel,
		 long gridPosition[NDIM], 
		 int nThreads=glb::num_omp_threads)
  {        
    typedef typename AMR::Voxel Voxel;
    //unsigned long ompDelta=std::max(1UL,count/nThreads/10);

#pragma omp parallel for num_threads(nThreads) //schedule(dynamic,ompDelta)
    for (long i=0;i<count;++i)
      {
	const Voxel *v=voxels[i];
	addVoxel<AMR>(v->getIndex(),v->getLevel(),v->data,
		      gridLevel,gridPosition);
      }    
  }

  template <class AMR>
  void addVoxels(const std::vector<typename AMR::Voxel *> &voxels, 
		 long gridLevel,
		 long gridPosition[NDIM], 
		 int nThreads=glb::num_omp_threads)
  {
    addVoxels<AMR>(&voxels[0],voxels.size(),gridLevel,gridPosition,nThreads);
  }

  template <class AMR>
  void addAmrGrid(AMR *amr, int nThreads=glb::num_omp_threads)
  {
    /*
    long level = AMR::MAX_LEVEL;
    long levelResolution = 1L<<level;
    long gridLevel=-1;

    long pos[NDIM];
    long parentResolution[NDIM];
    long resolution[NDIM];

    long newPos[NDIM];
    long newResolution[NDIM];

    if (params.haveParentGrid)
      {
	std::copy(params.position,params.position+NDIM,pos);
	std::copy(params.parentResolution,params.parentResolution+NDIM,parentResolution);
      }
    else
      {
	std::fill(pos,pos+NDIM,0);
	std::copy(params.resolution,params.resolution+NDIM,parentResolution);
      }

    std::copy(params.resolution, params.resolution+NDIM,resolution);
      
    for (int i=0;i<NDIM;++i)
      {
	
	long factor = levelResolution/parentResolution[i];
	if (factor*parentResolution[i] != levelResolution)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("This function is only implemented for cubic grids with sizes a power of 2.\n");
	    exit(-1);
	  }

	int p;
	for (p=0;(1<<p)<parentResolution[i];++p) {}
	if (gridLevel<0) gridLevel=p;
	else if (gridLevel!=p)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("This function is only implemented for cubic grids with size a power of 2\n");
	    exit(-1);
	  }

	newPos[i]=pos[i]*factor;
	newResolution[i]=resolution[i]*factor;
      }

    level -= AMR::ROOT_LEVEL;
    gridLevel-= AMR::ROOT_LEVEL;
    */
     
    long gridLevel;
    long pos[NDIM];
    long resolution[NDIM];
    
    long newPos[NDIM];
    long newResolution[NDIM];
    long newLevel;  
    
    addAmrGrid_getParams<AMR>(pos,newPos,resolution,newResolution,gridLevel,newLevel); 

    typedef typename AMR::Voxel Voxel;
    std::vector<Voxel*> voxels;
    amr->getLeavesGridOverlap(newPos,newResolution,newLevel,
			      std::back_inserter(voxels),true);	

    addVoxels<AMR>(voxels,gridLevel,pos,nThreads);

    // glb::console->print<LOG_STD_ALL>("OVERLAP: %ld cells overlap (out of %ld)\n",voxels.size(),amr->getNLeaves());
  }  

  template <class AMR>
  void addRemoteAmrGrid(const std::vector<typename AMR::Voxel::MpiExchangeData> &voxels,
			int nThreads=glb::num_omp_threads)
  {     
    addRemoteAmrGrid<AMR>(&voxels[0],voxels.size(),nThreads);
  }

  template <class AMR>
  void addRemoteAmrGrid(const typename AMR::Voxel::MpiExchangeData *voxels, 
			unsigned long count, 
			int nThreads=glb::num_omp_threads)
  {
    long gridLevel;
    long pos[NDIM];
    long resolution[NDIM];
    
    long newPos[NDIM];
    long newResolution[NDIM];
    long newLevel;  
    
    addAmrGrid_getParams<AMR>(pos,newPos,resolution,newResolution,gridLevel,newLevel); 
    addRemoteVoxels<AMR>(voxels,count,gridLevel,pos,nThreads);
  }  
  /*
  template <class InputIterator, class OutputIterator>
  void interpolateMany_CIC(InputIterator coords_start, InputIterator coords_stop,
			   OutputIterator result, int fieldIndex) const
  {
    const long delta = (fieldIndex<0)?nFields:1;
    double res[nFields];

    for (InputIterator it=coords_start; it!=coords_stop; ++it)
      {
	interpolate_CIC(*it,res,fieldIndex);
	for (int i=0;i<delta;++i) {*result=res[i];++result;}
      }
  }
  */

  template <int ORDER, class CT, bool CellCenteredValues=true> 
  void interpolate(CT *coords, double *result, int fieldIndex=-1) const
  {
    internal::KernelDispacherT<ORDER>::interpolate(*this,coords,result,fieldIndex);
  }

  template <class CT, bool CellCenteredValues=true> 
  void interpolate(int order, CT *coords, double *result, int fieldIndex=-1) const
  {
    switch (order)
      {
      case 0:
	internal::KernelDispacherT<0>::interpolate(*this,coords,result,fieldIndex);break;
      case 1:
	internal::KernelDispacherT<1>::interpolate(*this,coords,result,fieldIndex);break;
      case 2:
	internal::KernelDispacherT<2>::interpolate(*this,coords,result,fieldIndex);break;
      default:
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Order %d interpolation not implemented yet.\n",order);
	exit(-1);
	//do_interpolate<CellCenteredValues>(coords,result,fieldIndex,gridKernel_CIC);
      }
  }

  template <int ORDER, class CT, bool CellCenteredValues=true> 
  void interpolateGradient(CT *coords, double *result, int fieldIndex=-1) const
  {
    internal::KernelDispacherT<ORDER,1>::interpolate(*this,coords,result,fieldIndex);
  }

  template <class CT, bool CellCenteredValues=true> 
  void interpolateGradient(int order, CT *coords, double *result, int fieldIndex=-1) const
  {
    switch (order)
      {
      case 0:
	//internal::KernelDispacherT<0>::interpolate(*this,coords,result,fieldIndex);break;
      case 1:
	internal::KernelDispacherT<1,1>::interpolate(*this,coords,result,fieldIndex);break;
      case 2:
	internal::KernelDispacherT<2,1>::interpolate(*this,coords,result,fieldIndex);break;
      default:
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Order %d gradient interpolation not implemented yet.\n",order);
	exit(-1);
	//do_interpolate<CellCenteredValues>(coords,result,fieldIndex,gridKernel_CIC);
      }
  }

  /*
  // Checked only if value is defined at cells and for PBC!!!!!
  // fieldIndex<0 => all fields !
  template <class CT, bool CellCenteredValues=true> 
  void interpolate_TSC(CT *coords, double *result, int fieldIndex) const
  {    
    Data *ptr;
    long pos[NDIM];
    
    if (fieldIndex>=0)
      {
	Data* ptr = arr;
	
	if (IS_INTERLEAVED) ptr+=fieldIndex; 
	else ptr+=fieldIndex*getNValues();

	ptr += gridKernel_TSC.setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);

	gridKernel_TSC.apply<IS_PERIODIC,gridKernel::APPLY>(ptr,result,pos);
      }
    else
      {
	size_t delta= gridKernel_TSC.template setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);
	
	Data* ptr = arr+delta;
	for (int i=0;i<nFields;++i)
	  {	 
	    gridKernel_TSC.apply<IS_PERIODIC,gridKernel::APPLY>(ptr,result+i,pos);

	    if (IS_INTERLEAVED) ++ptr;
	    else ptr+=getNValues();	      
	  }
	//if (debug) printf("@(%e %e) => (%ld %ld) : %e %e\n",coords[0],coords[1],pos[0],pos[1],result[0],result[1]);
      }
  }
  
  // Checked only if value is defined at cells and for PBC!!!!!
  // fieldIndex<0 => all fields !
  template <class CT, bool CellCenteredValues=true> 
  void interpolate_CIC(CT *coords, double *result, int fieldIndex) const
  {
    
    Data *ptr;
    long pos[NDIM];
    
    if (fieldIndex>=0)
      {
	Data* ptr = arr;
	
	if (IS_INTERLEAVED) ptr+=fieldIndex; 
	else ptr+=fieldIndex*getNValues();

	ptr += gridKernel_CIC.setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);

	gridKernel_CIC.apply<IS_PERIODIC,gridKernel::APPLY>(ptr,result,pos);
      }
    else
      {
	size_t delta= gridKernel_CIC.template setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);

	Data* ptr = arr+delta;
	for (int i=0;i<nFields;++i)
	  {	 
	    gridKernel_CIC.apply<IS_PERIODIC,gridKernel::APPLY>(ptr,result+i,pos);

	    if (IS_INTERLEAVED) ++ptr;
	    else ptr+=getNValues();	      
	  }
      }

  
  }
  */
  
  template <class F>
  double quadrature(F &f, int field=0,
		    int nThreads=glb::num_omp_threads)
  {
    typedef GaussQuadratureT<NDIM> GQ;
    double result=0;
    
#pragma omp parallel num_threads(nThreads) reduction(+:result)    
      {
	int th=omp_get_thread_num();
	// We need a different copy for each thread as we do not know
	// if the functors are threadsafe !
	internal::lrg::QFunctorAdapter<F,NDIM> qFA(f,field);	

	oneFieldIterator it(begin(),field);
	const oneFieldIterator it_end(end(),field);
	it+=th;	
	for (;it!=it_end;it+=nThreads)
	  {
	    qFA.set(it);
	    result += GQ::compute(qFA);
	  }
      }

    return result;  
  }

  /*
    // NOT IMPLEMENTED YET
  template <class FL>
  void quadratureFromList(FL &functorList, double *result, int field=0,
			  int nThreads=glb::num_omp_threads)
  {
    typedef GaussQuadratureT<NDIM> GQ;
    static const int stride = 16; // To avoid false sharing (that's too large, for safety)
    static const int delta = (FL::SIZE+stride);
    std::vector<double> tmpResult(delta*nThreads,0);

#pragma omp parallel num_threads(nThreads) reduction(+:result)    
    {
      auto adaptedFunctorList = 
	internal::lrg::adaptFunctorList<NDIM>(functorList,field);
      int th=omp_get_thread_num();	
      double *res=&tmpResult[delta*th];

      oneFieldIterator it(begin(),field);
      const oneFieldIterator it_end(end(),field);
      it+=th;	
      for (;it!=it_end;it+=nThreads)
	{
	  double tmp[FL::SIZE];
	  internal::lrg::setAdaptedFunctorList(adaptedFunctorList,it);
	  GQ::computeList(adaptedFunctorList,tmp);
	  for (int i=0;i<FL::SIZE;++i) res[i]+=tmp[i];
	}
    }

    for (int i=1;i<nThreads;++i)
      for (int j=0;j<FL::SIZE;++j)
	{
	  tmpResult[j] += tmpResult[i*delta+j];
	}
      
    std::copy_n(&tmpResult[0],FL::SIZE,result);
  }
  */

  // Checked only if value is defined at cells and for PBC!!!!!
  // fieldIndex<0 => all fields !
  template <class K, class CT, bool CellCenteredValues=true> 
  void applyKernel(K& kernel, CT *coords, double *result, int fieldIndex=-1) const
  {    
    //Data *ptr;
    long pos[NDIM];
    
    if (fieldIndex>=0)
      {
	const Data* ptr = arr;
	
	if (IS_INTERLEAVED) ptr+=fieldIndex; 
	else ptr+=fieldIndex*getNValues();

	ptr += kernel.template setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);

	kernel.template apply<IS_PERIODIC,gridKernel::APPLY>
	  (ptr,result,pos);
      }
    else
      {
	size_t delta= kernel.template setCoefs<IS_PERIODIC,CellCenteredValues>
	  (coords,params.x0,boxSizeInv,pos);
	
	const Data* ptr = arr+delta;
	for (int i=0;i<nFields;++i)
	  {	 
	    kernel.template apply<IS_PERIODIC,gridKernel::APPLY>
	      (ptr,result+i*K::OUT_COUNT,pos);

	    if (IS_INTERLEAVED) ++ptr;
	    else ptr+=getNValues();	      
	  }
	//if (debug) printf("@(%e %e) => (%ld %ld) : %e %e\n",coords[0],coords[1],pos[0],pos[1],result[0],result[1]);
      }
  }
   
  void applySymmetry(RegularGridSymmetryE symmetry,
		    int nThreads=glb::num_omp_threads)
  {
    //unsigned long s=static_cast<unsigned long>(symmetry);
    internal::lrg::SymmetryFunctorInterface<MyType> *functor = nullptr;

    switch (symmetry)
      {
      case RGSV::NONE:
	return;
	break;

      case RGSV::PLANAR0:
	functor = new internal::lrg::SymmetryFunctorSinglePlane<MyType,0>(this);
	break;
      case RGSV::PLANAR1:
	functor = new internal::lrg::SymmetryFunctorSinglePlane<MyType,1>(this);
	break;
      case RGSV::PLANAR2:
	functor = new internal::lrg::SymmetryFunctorSinglePlane<MyType,2>(this);
	break;
      case RGSV::PLANAR01:
	functor = new internal::lrg::SymmetryFunctorDualPlane<MyType,0,1>(this);
	break;
      case RGSV::PLANAR02:
	functor = new internal::lrg::SymmetryFunctorDualPlane<MyType,0,2>(this);
	break;
      case RGSV::PLANAR12:
	functor = new internal::lrg::SymmetryFunctorDualPlane<MyType,1,2>(this);
	break;
      case RGSV::PLANAR012:
	functor = new internal::lrg::SymmetryFunctorTriplePlane<MyType>(this);
	break;  

      case RGSV::CONSTANT0:
	functor = new internal::lrg::SymmetryFunctorConstant<MyType,0>(this);
	break;
      case RGSV::CONSTANT1:
	functor = new internal::lrg::SymmetryFunctorConstant<MyType,1>(this);
	break;
      case RGSV::CONSTANT2:
	functor = new internal::lrg::SymmetryFunctorConstant<MyType,2>(this);
	break;

      case RGSV::CENTRAL:
	functor = new internal::lrg::SymmetryFunctorCentral<MyType>(this);
	break;  
      default:	
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("Unimplemented symmetry type : %s(=%ld) \n",
	   RegularGridSymmetrySelect().getString(symmetry).c_str(),
	   static_cast<long>(symmetry));
      }

#pragma omp parallel num_threads(nThreads)
    {
      int th=omp_get_thread_num();

      inOrderIterator it(begin());
      const inOrderIterator it_end(end());
      it+=th;	  

      int N=functor->getMaxCount();
      Data *dest[N]; 
      for (;it!=it_end;it+=nThreads)
	{
	  int n=functor->getSymmetries(it,dest);
	  if (n>0)
	    {
	      for (int i=0;i<n;++i) (*it)+=*dest[i];
	      (*it)/=(n+1);
	      for (int i=0;i<n;++i) *dest[i]=(*it);
	    }
	}
      }
    delete functor;
  }

  
private:  
  template <class AMR, class T>
  void addAmrGrid_getParams(T pos[NDIM], T newPos[NDIM], 
			    T resolution[NDIM], T newResolution[NDIM],
			    T&level,T &newLevel)
  {
    newLevel = AMR::MAX_LEVEL;
    long levelResolution = 1L<<newLevel;
    level=-1;
   
    long parentResolution[NDIM];

    if (params.haveParentGrid)
      {
	std::copy(params.position, params.position+NDIM, pos);
	std::copy(params.parentResolution, params.parentResolution+NDIM, parentResolution);
      }
    else
      {
	std::fill(pos,pos+NDIM,0);
	std::copy(params.resolution, params.resolution+NDIM, parentResolution);
      }

    std::copy(params.resolution, params.resolution+NDIM, resolution);
      
    for (int i=0;i<NDIM;++i)
      {	
	long factor = levelResolution/parentResolution[i];
	if (factor*parentResolution[i] != levelResolution)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("This function is only implemented for cubic grids with sizes a power of 2.\n");
	    exit(-1);
	  }

	int p;
	for (p=0;(1<<p)<parentResolution[i];++p) {}
	if (level<0) level=p;
	else if (level!=p)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("This function is only implemented for cubic grids with size a power of 2\n");
	    exit(-1);
	  }

	newPos[i]=pos[i]*factor;
	newResolution[i]=resolution[i]*factor;
      }

    newLevel -= AMR::ROOT_LEVEL;
    level-= AMR::ROOT_LEVEL;    
  }

  /*
    NDF::NDfield *toNDfield()
    {
    long i;
    const Params &gp=getParams();
    int useSpecies=(gp.nFields>1)?1:0;
    int dims[P_NDIM+U_NDIM+useSpecies];
    double x0[P_NDIM+U_NDIM+useSpecies];
    double delta[P_NDIM+U_NDIM+useSpecies];
    int ndims;
    long stride[P_NDIM+U_NDIM+1+useSpecies];
    int index;
    value_type *d;
    char comment[80];
    //gridT *grid=bT::gh->getGrid();
    

    stride[0]=1;    
    if (useSpecies)
    {
    dims[0]=gp.nFields;
    x0[0]=0;
    delta[0]=gp.nFields;
    stride[1]=stride[0]*dims[0];
    }

    long dp=useSpecies;
    for (i=0;i<NDIM;i++) {
    dims[i+dp]=getValCoord(i).size();
    x0[i+dp]=getValCoord(i).front();
    delta[i+dp]=getValCoord(i).back()-getValCoord(i).front();
    stride[i+1+dp]=getValStride(i+1);
    }
    
    strcpy(comment,"");
    for (i=0;i<P_NDIM;i++) 
    {
    char tmp[256];
    sprintf(tmp,"%d %d ",gparams.scale[i],(int)valLocation_P);
    strcat(comment,tmp);
    }
    for (i=0;i<U_NDIM;i++) 
    {
    char tmp[256];
    sprintf(tmp,"%d %d ",gp.Uscale[i],(int)valLocation_U);
    strcat(comment,tmp);
    }
    
   
    d = arr;
    ndims=P_NDIM+U_NDIM+useSpecies;
 
    return NDF::Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,d,comment);
    }
    
    void writeToNDfield(const std::string &fname_)
    {
    char fname[255];
    sprintf(fname,"%s.ND",fname_.c_str());
    
    NDF::NDfield *f=toNDfield();
    NDF::Save_NDfield(f,fname);
    free(f);
    }
  */
 
protected:
  /*
  template <class AMR, class G>
  class AmrGridVisitorT
  {
    typedef typename AMR::Voxel Voxel;  
    static const in NDIM = G::NDIM;

  public:
    AmrGridVisitorT(AMR *amr_, G *grid):amr(amr_)
    {
      grid->getBoundingBox(bbox);
    }

    bool visit(Voxel *voxel) const
    {
      if (voxel->isLeaf())
	{
	  

	  return false;
	}
      else
	{
	  index2CornerCoordsAndOpp
	}

      return true;
    }    

    static void visited(const Voxel *voxel) 
    {
    
    }
    
  private:
    
    const AMR *amr; 
    double bbox[NDIM][2];
  };
  */
};

/** \}*/
#include "../internal/namespace.footer"
#endif
