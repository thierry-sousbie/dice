#ifndef __REGULAR_GRID_HXX__
#define __REGULAR_GRID_HXX__

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <vector>

#include "../dice_globals.hxx"

#include "./valLocationType.hxx"
#include "./localRegularGrid.hxx"
#include "./gridTopology.hxx"
#include "./regularGridSlicer_simple.hxx"
#include "./gridKernels.hxx"

#include "../geometry/boundaryType.hxx"
#include "../geometry/geometricProperties.hxx"

#include "../IO/vtkPRectilinearGridWriter.hxx"

#include "../tools/helpers/helpers.hxx"
#include "../tools/IO/myIO.hxx"
#include "../tools/IO/paramsParser.hxx"
#include "../tools/MPI/mpiCommunication.hxx"

/**
 * @file 
 * @brief A regular grid shared among MPI processes
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

/**
 * \class RegularGridT
 * \brief A regular grid shared among MPI processes that manages local grid chunks of
 * type LocalRegularGridT.
 * \tparam ND The number of dimensions
 * \tparam DT The type of data to store at grid sites
 * \tparam BT The boudary conditions type (see BoundaryType)
 * \tparam FML defines the memory layout of the different fields in the grid ( see 
 * regularGridFieldLayout )
 * \note Data is stored in column-major order (i.e. Fortran order) : the first dimension
 * index changes first when traversing the array linearly in memory.
 */
template <long ND, 
	  typename DT = double,
	  long BT = BoundaryType::NONE,
	  long FML = regularGridFieldLayout::CONSECUTIVE>
class RegularGridT {
public:

  typedef RegularGridT<ND,DT,BT,FML> MyType;
  typedef LocalRegularGridT<ND,DT,BT,FML> LocalGrid;
 
  static const long BOUNDARY_TYPE = BT;
  static const int  NDIM = ND;

  static const int  IS_INTERLEAVED = (FML==regularGridFieldLayout::INTERLEAVED);
  static const int  IS_PERIODIC = (BT==BoundaryType::PERIODIC);
  static const int  LAYOUT = FML;

  typedef DT value_type;
  typedef DT Data;
  
  static std::string classHeader() {return "regular_grid";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  typedef RegularGridSlicerSimpleT<MyType> DefaultSlicer; /**< a default grid slicer */

  typedef typename LocalGrid::Scale Scale;
  typedef typename LocalGrid::ScaleType ScaleType;
  typedef typename LocalGrid::ScaleTypeV ScaleTypeV;
  typedef typename LocalGrid::Params Params;
  typedef typename LocalGrid::GridNav GridNav;
  typedef typename LocalGrid::Direction Direction;

  typedef typename LocalGrid::iterator local_iterator;
  typedef typename LocalGrid::fieldIterator local_fieldIterator;
  typedef typename LocalGrid::oneFieldIterator local_oneFieldIterator;
  typedef typename LocalGrid::inOrderIterator local_inOrderIterator;
  
  //typedef typename LocalGrid::ValLocationType ValLocationType;
  //typedef typename LocalGrid::ValLocationTypeV ValLocationTypeV;

  typedef GridTopologyT<NDIM> GridTopology;

protected:
  std::vector< std::pair<long,Direction> > sendDir;  
  std::vector< std::pair<local_fieldIterator,local_fieldIterator> > sendIterators;
  std::vector< std::vector<Data> > sendBuffer;

  std::vector< std::pair<long,Direction> > receiveDir;  
  std::vector< std::pair<local_fieldIterator,local_fieldIterator> > receiveIterators;
  std::vector< std::vector<Data> > receiveBuffer;
  
  // setup grid boundaries synchronization ...
  // Used only when local grids are overlapping, used to work but not checked
  // for a long time, may have to rewrite this if usefull ...
  void initSync()
  {
    long i,j,k;

    std::vector< std::vector<Direction> > marginDir(NDIM);
    std::vector< long > buffer((2+NDIM)*pow(3.0,double(NDIM)));
 
    for (i=0;i<NDIM;i++) 
      {
	marginDir[i].push_back(GridNav::dir(i,0));
	if (grid->getLowMargin(i)>0)  marginDir[i].push_back(GridNav::dir(i,-1));
	if (grid->getHighMargin(i)>0) marginDir[i].push_back(GridNav::dir(i,+1));
      }
    
    std::vector<int> w(NDIM,0);
    while (true)
      {
	Direction direction=0;
	int delta=1;
	for (i=0;i<NDIM;i++) 
	  {
	    direction |= marginDir[i][w[i]];
	    w[i]+=delta;
	    if (w[i]>=marginDir[i].size()) {w[i]=0;delta=1;}
	    else delta=0;
	  }

	if (direction) 
	  {	   
	    long nei = getNeighborIndex(direction);
	    receiveDir.push_back(std::make_pair(nei,direction));
	   
	    local_iterator bit=grid->margin_begin(receiveDir.back().second);
	    local_iterator eit=grid->margin_end(receiveDir.back().second);
	    
	    receiveIterators.push_back
	      (std::make_pair(local_fieldIterator(bit),local_fieldIterator(eit)));
	  }

	if (delta) break;
      }
   
    receiveBuffer.resize(receiveDir.size());
    for (i=0;i<receiveDir.size();i++)
      receiveBuffer[i].resize(grid->getNFields()*grid->margin_size(receiveDir[i].second));

    for (i=0;i<mpiCom->size();i++)
      {
	if (i==mpiCom->rank())
	  {
	    buffer[0]=receiveDir.size();
	    for (j=0;j<receiveDir.size();j++)
	      {
		buffer[1+(2+NDIM)*j]=(long)receiveDir[j].first;
		buffer[1+(2+NDIM)*j+1]=(long)GridNav::reverse(receiveDir[j].second);
		local_iterator it=grid->margin_begin(receiveDir[j].second);
		for (k=0;k<NDIM;k++) buffer[1+(2+NDIM)*j+2+k]=local_iterator::dim(it,k);
	      }
	  }

	mpiCom->Bcast(buffer,i);
	//MPI_Bcast(&buffer[0],buffer.size(),MPI_LONG,i,MPI_COMM_WORLD);

	int sendInfo[NDIM];
	for (j=0;j<buffer[0];j++)
	  {
	    if (buffer[1+(2+NDIM)*j] == mpiCom->rank())
	      {
		sendDir.push_back(std::make_pair(i,buffer[1+(2+NDIM)*j+1]));
		for (k=0;k<NDIM;k++) sendInfo[k]=buffer[1+(2+NDIM)*j+2+k];	
		//printf("send : %d %d %d %ld\n",sendInfoP[0],sendInfoU[0],sendInfoU[1],sendDir.back().second);	
		local_iterator bit=grid->innerMargin_begin(sendInfo,sendDir.back().second);
		local_iterator eit=grid->innerMargin_end(sendInfo,sendDir.back().second);
	
		sendIterators.push_back
		  (std::make_pair(local_fieldIterator(bit),local_fieldIterator(eit)));
	      }
	  }
      }
   
    sendBuffer.resize(sendDir.size());
    for (i=0;i<sendDir.size();i++)
      sendBuffer[i].resize(grid->getNFields()*
			   local_iterator::boxSize(sendIterators[i].first));      
    
    if (glb::debug)
      {
	char tmpS[10000];
	char tmpR[10000];

	strcpy(tmpS,"");
	strcpy(tmpR,"");

	for (i=0;i<sendDir.size();i++) 
	  sprintf(tmpS,"%s (%ld,%ld:%ld)",tmpS,sendDir[i].first,(long)sendDir[i].second,
		  local_iterator::boxSize(sendIterators[i].first));
	for (i=0;i<receiveDir.size();i++) 
	  sprintf(tmpR,"%s (%ld,%ld:%ld)",tmpR,receiveDir[i].first,
		  (long)receiveDir[i].second,
		  local_iterator::boxSize(receiveIterators[i].first));

	glb::console->print<LOG_DEBUG>
	  ("  (%d:) I send data to : %s.\n  (%d:) I receive data from: %s\n",
	   mpiCom->rank(),tmpS,mpiCom->rank(),tmpR);
      }  
  }

  void init()
  {
    initSync();

    for (int i=0;i<NDIM;++i)
      {
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

	if (grid->getParams().x0[i]<=params.x0[i]) 
	  onLowBoundary[i]=1;
	else
	  onLowBoundary[i]=0;

	if (grid->getParams().x0[i]+grid->getParams().delta[i] >=
	    params.x0[i] + params.delta[i]) 
	  onHighBoundary[i]=1;
	else
	  onHighBoundary[i]=0;

	onBoundary[i]=(onHighBoundary[i]||onLowBoundary[i]);
      }

    initialized=true;    
  }

  bool isInitialized() const {return initialized;}

public:

  /** \brief Save the grid to a VTK PRectilinearGrid file
   *  \param globalFName The name of the global file (no extension)
   *  \param format The format of the name for local files (no extension). This should 
   *  contain a '%d' that will be replaced by the process rank
   */
  void toVtk(const char *globalFName, const char *format)
  {
    IO::VtkPRectilinearGridWriterT<MyType> vtkWriter(this,mpiCom,globalFName,format);
    vtkWriter.write();
  }

  /** \brief initialize a finite difference / interpolation kernel (see gridKernel)
   */
  template <class K>
  void intializeKernel(K &kernel) const
  {
    grid->initialize(kernel);
  }
  
  /** \brief synchronize the boundaries of the local grids if any
   */
  void synchronize()
  {
    MPI_Request requestS[sendDir.size()];
    MPI_Request requestR[receiveDir.size()];
    long i;
    
#pragma omp parallel for 
    for (i=0;i<sendDir.size();i++)
      std::copy(sendIterators[i].first,sendIterators[i].second,sendBuffer[i].begin());
     
    for (i=0;i<receiveDir.size();i++)
      mpiCom->Irecv(&receiveBuffer[i][0],receiveBuffer[i].size(),receiveDir[i].first,&requestR[i],0);

    for (i=0;i<sendDir.size();i++)
      mpiCom->Isend(&sendBuffer[i][0],sendBuffer[i].size(),sendDir[i].first,&requestS[i],0);
   
    mpiCom->Waitall(requestS);
    mpiCom->Waitall(requestR);
   
#pragma omp parallel for   
    for (i=0;i<receiveDir.size();i++)
      std::copy(receiveBuffer[i].begin(),receiveBuffer[i].end(),receiveIterators[i].first);
   }
  
  //! constructor
  explicit RegularGridT(const char *dataName_="value"):
    grid(NULL),
    dataName(dataName_),
    initialized(false)
  {
    numThreads=glb::num_omp_threads;
  }
  
  //! destructor
  ~RegularGridT()
  {
    if (grid!=NULL) delete grid;
  }  

  void clone(MyType &cloned, bool cloneExtraElements=true) const
  { 
    LocalGrid *tmpGrid;

    if (!initialized) return;

    if (!cloned.initialized)
      tmpGrid=new LocalGrid();
    else
      tmpGrid=cloned.grid;    
    
    cloned = *this;
    cloned.grid = tmpGrid;
    grid->clone(*cloned.grid,cloneExtraElements);
  }

  template <class Source>
  long copyRawBuffer(Source &source)
  {
    return grid->copyRawBuffer(*source.grid);
  }

  RegularGridT( const MyType& other ):
    grid(NULL),
    dataName("copyConstructed"),
    initialized(false)
  {
    other.clone(*this);
  }
private:
  // Prevent public copy but keep default copy for clone implementation
  RegularGridT& operator=(const MyType&) = default;

public:
  /*
  void setNumThreads(int nThreads)
  {
    numThreads=nThreads;
  }
  */

  /*
  template <class Slicer = DefaultSlicer>
  void initialize(const Params &gp,
		  MpiCommunication *com=glb::mpiComWorld,
		  int nThreads = glb::num_omp_threads)
  {
    Slicer slicer;
    initialize(gp,slicer,com,nThreads);
  }
  */

  /** \brief initialize the grid according according to a slicer objet.
   *  See DefaultSlicer or RegularGridSlicerSimpleT for instance.
   */
  template <class Slicer>
  void initialize(Slicer &slicer, value_type* localGridArr=NULL)
  {
    if (initialized)
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>("The grid was initialized more than once !\n");
      }
 
    params=slicer.getGlobalGridParams();    
    mpiCom=slicer.getMpiCom();
    nGrids=slicer.getNChunks();   
    numThreads=slicer.getNThreads();

    //slicer.slice(params, mpiCom->rank(), mpiCom, nThreads);    
    slicer.slice(mpiCom->rank());
    for (int i=0;i<NDIM;++i)
      {
	sliceCount[i]=slicer.getSliceCount(i);
	slicePos[i]=slicer.getSlicePos(i);
	//localGridPos[i]=slicer.getLocalGridPos(i);
      }  
   
    if (grid!=NULL) delete grid;
    grid=new LocalGrid();
    grid->setName(dataName.c_str());
    grid->initialize(slicer.getLocalGridParams(),localGridArr);

    remoteData.resize(mpiCom->size());    
    topology.reset(mpiCom->size());
    for (int i=0;i<mpiCom->size();++i)
      {
	//slicer.slice(params, i, mpiCom, nThreads);
	slicer.slice(i);
	topology.setNeighbors(i,slicer);
	remoteData[i].set(slicer);
      }
    
    init();   
  } 
  
  template <class R>
  void initialize(R* reader,
		  MpiCommunication *com=glb::mpiComWorld,
		  int nThreads=glb::num_omp_threads)
  {
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);

    Params gp;    
    gp.read(reader);

    params=gp;
    mpiCom=com;
    numThreads=nThreads;

    int tmp;
    reader->read(&tmp);
    if (tmp!=NDIM)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Reading from binary file : Dimensions differ.\n");
      }
    reader->read(&tmp);
    if (tmp!=sizeof(Data))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Reading from binary file : Data type differ.\n");
      }
    reader->read(&nGrids);
    if (nGrids!=mpiCom->size())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Reading from binary file : mpi communication size differ.\n");
      }
    reader->read(&tmp);
    if (tmp!=mpiCom->rank())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Reading from binary file : mpi rank differ.\n");
      }
    reader->read(dataName);
    reader->read(sliceCount,NDIM);
    reader->read(slicePos,NDIM);
    //reader->read(localGridPos,NDIM);
    
    remoteData.resize(mpiCom->size()); 
    for (unsigned long i=0;i<remoteData.size();++i) 
      remoteData[i].read(reader);
   
    topology.read(reader);
  
    grid=new LocalGrid();
    grid->initialize(reader);

    init();    
  }

  template <class W>
  void write(W *writer)
  {
    writer->writeHeader(classHeader(),classVersion());
    params.write(writer);
    
    int d=NDIM;
    int sz=sizeof(Data);
    int rk=mpiCom->rank();

    writer->write(&d);
    writer->write(&sz);
    writer->write(&nGrids);
    writer->write(&rk);
    writer->write(dataName);

    writer->write(sliceCount,NDIM);
    writer->write(slicePos,NDIM);
    //writer->write(localGridPos,NDIM);

    for (unsigned long i=0;i<remoteData.size();++i)
      remoteData[i].write(writer);
    //remoteData[i].params.write(writer);
      
    topology.write(writer);
    //writer->write(neighborIndex,NDIM*2);    

    grid->write(writer);
  }

  std::string getName() const
  {
    return dataName;
  }

  void setName(const char *name)
  {
    dataName=name;
    if (grid!=NULL)
      grid->setName(name);
  }

  template <class K>
  void initializeInterpolationKernel(K &kernel, bool setSpacing=true) const
  {
    grid->initializeInterpolationKernel(kernel,setSpacing);
  }

  template <class K, class CT, bool CellCenteredValues=true> 
  void applyKernel(K& kernel, CT *coords, double *result, int fieldIndex=-1) const
  {
    grid->template applyKernel<K,CT,CellCenteredValues>(kernel,coords,result,fieldIndex);
  }
 
  Data *getDataPtr() {return grid->getDataPtr();}
  value_type *getDataPtr(long i) {return grid->getDataPtr(i);}
  value_type *getDataPtr(long i, int field) {return grid->getDataPtr(i,field);}
  template <typename T>
  value_type *getDataPtr(T* w) {return grid->getDataPtr(w);}
  template <typename T>
  value_type *getDataPtr(T* w, int field) {return grid->getDataPtr(w,field);}

  LocalGrid *getLocalGrid() {return grid;}

  const Params &getGlobalGridParams() const {return params;}
  const Params &getLocalGridParams() const {return grid->getParams();}  
  const Params &getRemoteGridParams(int index) const {return remoteData[index].params;}  

  const std::vector<double>& getCellCoord(long i) const 
  {return cellCoord[i];}
  const std::vector<double>& getVertexCoord(long i) const 
  {return vertexCoord[i];}
  const std::vector<double>& getValueCoord(long i) const 
  {return valueCoord[i];}  

  double getConstantVolumeElement() const
  {return grid->getConstantVolumeElement();}
  
  double getSize(int dim) const {return params.delta[dim];}
  double getLocalSize(int dim) const {return grid->getSize(dim);}

  double getOrigin(int dim) const {return params.x0[dim];}
  double getLocalOrigin(int dim) const {return grid->getOrigin(dim);}

  //long getSize(int dim) const {return valueCoord[dim].size();}
  //long getLocalSize(int dim) const {return grid->getSize(dim);}
  long getResolution(int dim) const {return params.resolution[dim];}
  const long *getResolution() const {return params.resolution;}
  long getLocalResolution(int dim) const {return grid->getResolution(dim);}
  const long *getLocalResolution() const {return grid->getResolution();}
  long getLocalValueStride(int dim) const {return grid->getValueStride(dim);}
  const long *getLocalValueStride() const {return grid->getLocalValueStride();}
  long getNFields() const {return grid->getNFields();}
  long getLocalNValues() const {return grid->getNValues();}

  ScaleType getScaleType(int i) const {return params.scale[i];}
  ValLocationType getValLocation(int dim) const {return params.valLocation[dim];}

  bool onGlobalLowBoundary(int dim) const {return onLowBoundary[dim];}
  bool onGlobalHighBoundary(int dim) const {return onHighBoundary[dim];}
  bool onGlobalBoundary(int dim) const {return onBoundary[dim];}

  int getSliceCount(int dim) const {return sliceCount[dim];}
  int getSlicePos(int dim) const {return slicePos[dim];}
  int getLocalPosition(int dim) const {return grid->getParams().position[dim];}  
  int getLocalPosition(int gridIndex, int dim) const 
  {return remoteData[gridIndex].params.position[dim];}

  template <class CT>
  void getBoundingBox(CT bbox[NDIM][2]) const
  {
    for (int i=0;i<NDIM;++i)
      {
	bbox[i][0] = params.x0[i];
	bbox[i][1] = params.x0[i] + params.delta[i];
      }
  }

  void reinterpretNFields(int N=0)
  {
    grid->reinterpretNFields(N);
  }
  
  int getRemoteNeighborIndex(int index, int dim, int dir) const
  {
    return topology.getNeighbor(index,dim,dir);
  }

  int getNeighborIndex(int dim, int dir) const
  {
    return topology.getNeighbor(mpiCom->rank(),dim,dir);    
  }

  int getRemoteNeighborIndex(int index, Direction dir) const
  {
    return topology.getNeighbor(index,dir);
  }

  int getNeighborIndex(Direction dir) const
  {
    return topology.getNeighbor(mpiCom->rank(),dir);
  }

  /*
  int getNeighborIndex(int dim, int dir) const
  {
    return topology.getNeighbor(mpiCom->rank(),dim,dir);
    // if (dir==0) 
    //   return -1;
    // else 
    //   return neighborIndex[dim][(dir>0)];
  }

  int getNeighborIndex(Direction dir) const
  {
    for (int i=0;i<NDIM;++i)
      {
	if (dir == GridNav::dir(i,-1))
	  return neighborIndex[i][0];
	else if (dir == GridNav::dir(i,+1))
	  return neighborIndex[i][1];
      }
    return -1;
  }
  */

  /*
  template <class SubGrid>
  void getSubGrid(SubGrid &result,int which[SubGrid::NDIM], 
		  int nFields, bool global=true, long minStorage=0)
  {
    typename SubGrid::Params subP;

    if (global) 
      getGlobalGridParams().getSubGridParams(subP,which,minStorage);
    else
      grid->getParams().getSubGridParams(subP,which,minStorage);
    if (nFields>=1) subP.nFields=nFields;
    
    result.initialize(subP);
  }
  */

  /*
  template <class iterator>
  std::vector<double> gather_all(const iterator &tab_start, long dir, long addSpare=0)
  {
    int dirSize=getGlobalGridSize(dir);
    std::vector<double> result(dirSize+addSpare,0);
    long i;   
    
    std::vector<double> Total_in(getSliceCount(dir)+addSpare,0);
    std::vector<double> Total_out(getSliceCount(dir)+addSpare,0);
    
    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValueCoord(dir).size()-grid->getHighMargin(dir);       
    
    int pos=getLocalGridPos(dir);
  
    std::copy(tab_start+imin,tab_start+imax,&result[addSpare+pos+imin]);
    mpiCom->Allreduce_inplace(result,MPI_SUM);
   
    return result;
  }
  */

  /*
  template <class iterator>
  void gather(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> tmp=gather_all(tab_start,dir,addSpare);

    int pos=getLocalGridPos(dir);
    long size = grid->getValueCoord(dir).size();
    
    std::copy(&tmp[pos],&tmp[pos+size],tab_start);    
  }
  */
  /*
  template <class iterator>
  std::vector<double> accumulate_all(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> result=gather_all(tab_start,dir,addSpare);
    long i;
   
    for (i=1;i<result.size();i++) result[i]+=result[i-1];
      
    return result;
  }
  */
  /*
  template <class iterator>
  double accumulate(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> tmp=accumulate_all(tab_start,dir,addSpare);

    int pos=getLocalGridPos(dir);
    long size=grid->getValueCoord(dir).size();    

    std::copy(&tmp[pos],&tmp[pos+size],tab_start);
    return tmp.back();
  }
  */
  /*
  template <class iterator>
  double sum(const iterator &tab_start, long dir)
  {
    long i;   
    double sum=0;
    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValueCoord(dir).size()-grid->getHighMargin(dir);
    iterator it=tab_start+imin;
    for (i=imin;i<imax;i++,++it) sum+=(*it);
      
    mpiCom->Allreduce_inplace(&sum,1,MPI_SUM);
    
    return sum;
  }
  */
  /*
  double sum(double val)
  {
    double res=val;
    mpiCom->Allreduce_inplace(&res,(long)1,MPI_SUM);   
    return res;
  }
  */

  double getMax()
  {    
    double min,max;
    getLocalGrid()->getMinMax(min,max);
    return mpiCom->max(max);
  }
 
  template <class Functor>
  void visit(const Functor &f, int field=0, int nThreads=glb::num_omp_threads)
  {
    getLocalGrid()->visit(f,field,nThreads);
  }

  void subMultiply(double sub, double mult)
  {
    getLocalGrid()->subMultiply(sub,mult);
  }

  void erase()
  {
    grid->erase();
  }

  template <class AMR>
  void addAmrGrid(AMR *amr, int nThreads=glb::num_omp_threads, bool quiet=false)
  {
    typename TimerPool::Timer timer("dummy");

    typedef typename AMR::Voxel Voxel;
    typedef typename Voxel::MpiExchangeData MpiData;

    MpiDataType mpiDataType = MpiData::createMpiStructType();

    const long sz=mpiCom->size();
    const long rk=mpiCom->rank();

    //std::vector<Voxel*> voxels[sz];
    std::vector< std::vector<Voxel*> > voxels(sz);
    std::vector<int> sendCount(sz);

    //if (!quiet) glb::console->printFlush<LOG_STD>("Converting AMR to regular grid ... ");    
    // First take care of the local AMR
    sendCount[rk]=0;
    if (!quiet) glb::console->printFlush<LOG_STD>("    * Painting local chunk ... ");
    timer.start();
    grid->addAmrGrid(amr,nThreads);
    double elapsed = timer.stop();
    if (!quiet) glb::console->printFlush<LOG_STD>("done in %.3gs\n",elapsed);
    
    if (sz==1) return;

    if (!quiet) glb::console->printFlush<LOG_STD>
		  ("    * Preparing chunks for sending ... ");
    timer.start();
    // and then get the overlapping voxels of the (remote) AMR grids with each regular grid 
#pragma omp parallel for num_threads(nThreads) 
    for (int i=0;i<sz;++i)
      {		
	if (i!=rk)
	  {
	    long pos[NDIM];
	    long resolution[NDIM];
	    long level;	    

	    addAmrGrid_getParams<AMR,long>(i,pos,resolution,level);	    	    
	    amr->getLeavesGridOverlap(pos,resolution,level,
				      std::back_inserter(voxels[i]),true);	
	    sendCount[i]=(int)voxels[i].size();
	    
	    // glb::console->print<LOG_STD_ALL>
	    //   ("REMOTE OVERLAP(%d): %ld cells overlap (out of %ld)\n",
	    //    i,voxels[i].size(),amr->getNLeaves());	    
	  }
      }
    elapsed = timer.stop();
    if (!quiet) glb::console->printFlush<LOG_STD>("done in %.3gs\n",elapsed);

    if (!quiet) glb::console->printFlush<LOG_STD>("    * Exchanging everything ... ");
    timer.start();
    
    std::vector<int> receiveCount(sz);
    mpiCom->Alltoall(sendCount,receiveCount);    
    
    // USING Alltoallv
    std::vector<int> sendDisp(sz+1);
    sendDisp[0]=0;
    std::copy(sendCount.begin(),sendCount.end(),&sendDisp[1]);
    
    std::vector<int> receiveDisp(sz+1);
    receiveDisp[0]=0;
    std::copy(receiveCount.begin(),receiveCount.end(),&receiveDisp[1]);
    
    for (unsigned long i=1;i<(sz+1);++i)
      {
	sendDisp[i]+=sendDisp[i-1];
	receiveDisp[i]+=receiveDisp[i-1];
      }

    std::vector<MpiData> sendData(sendDisp.back());
    std::vector<MpiData> receiveData(receiveDisp.back());
    
    //#pragma omp parallel for num_threads(nThreads) 
    for (int i=0;i<sz;++i)
      {
	if (sendCount[i]>0)
	  std::copy(voxels[i].begin(),voxels[i].end(),&sendData[sendDisp[i]]);
      }
    
    MPI_Datatype dataType = mpiDataType.getType();
    mpiCom->Alltoallv(&sendData[0],&sendCount[0],&sendDisp[0],dataType,
		      &receiveData[0],&receiveCount[0],&receiveDisp[0],dataType);

    elapsed = timer.stop();
    if (!quiet) glb::console->printFlush<LOG_STD>
		  ("done in %.3gs\n",elapsed);
    if (!quiet) glb::console->printFlush<LOG_STD>
		  ("    * Painting remote chunks ... ");
    timer.start();
    for (int i=0;i<sz;++i)
      {
	// glb::console->printFlush<LOG_STD>
	//   ("Adding %d voxels from process %d.\n",receiveCount[i],i);
	if (receiveCount[i]>0)
	  grid->template addRemoteAmrGrid<AMR>
	    (&receiveData[receiveDisp[i]],receiveCount[i],nThreads);    
      }

    /* ------------------------------- */
    
    /*
    std::vector<MpiData> sendData[sz];
    std::vector<MpiData> receiveData[sz];
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> receiveRequests;

    sendRequests.reserve(sz);
    receiveRequests.reserve(sz);
    
    // Exchange the overlapping voxels 
    // FIXME: add OMP and protect critical regions ? (really usefull ?)
    //#pragma omp parallel for num_threads(nThreads) 
    for (int i=0;i<sz;++i)
      {
	if (sendCount[i]>0)
	  {
	    // glb::console->print<LOG_STD_ALL>("Process %d sending %d to %d\n",
	    //  				     rk,sendCount[i],i);
	    sendData[i].reserve(sendCount[i]);
	    sendData[i].assign(voxels[i].begin(),voxels[i].end());
	    long reqId=sendRequests.size();
	    sendRequests.resize(reqId+1);
	    mpiCom->IsendMpiType(&sendData[i][0],sendData[i].size(),
				 mpiDataType.getType(),i,
				 &sendRequests[reqId]);
	  }
	if (receiveCount[i]>0)
	  {
	    // glb::console->print<LOG_STD_ALL>("Process %d receiving %d from %d\n",
	    //  				     rk,receiveCount[i],i);
	    receiveData[i].resize(receiveCount[i]);
	    long reqId=receiveRequests.size();
	    receiveRequests.resize(reqId+1);
	    mpiCom->IrecvMpiType(&receiveData[i][0],receiveData[i].size(),
				 mpiDataType.getType(),i,
				 &receiveRequests[reqId]);
	  }
      }
    
    // Receive the overlapping voxels and process them. We cannot do that in parallel
    // because remote voxels may overlap ...    
    for (unsigned long i=0;i<receiveRequests.size();++i)
      {	
	MPI_Status status;
	int reqId;
	unsigned long count = receiveRequests.size()-i;

	mpiCom->Waitany(count,&receiveRequests[0],&reqId,&status);
	std::swap(receiveRequests[reqId],receiveRequests[count-1]);
	
	// process the request
	int reqRk = status.MPI_SOURCE;
	// glb::console->print<LOG_STD_ALL>("Process %d adding remote grid %d (%ld el.)\n",
	//  				 rk,index,receiveData[index].size());
	grid->template addRemoteAmrGrid<AMR>(receiveData[reqRk],nThreads);
      }        

    // Need to be sure everything was sent before we deallocate the buffers 
    mpiCom->Waitall(sendRequests);
    */
    elapsed = timer.stop();
    if (!quiet) glb::console->printFlush<LOG_STD>("done in %.2gs\n",elapsed);
    //if (!quiet) glb::console->printFlush<LOG_STD>("done.\n");
  }
  /*
  template <class DT>
  void gatherSubBoxC(DT pMin[NDIM], DT pMax[NDIM], 
				LocalGrid &localSubBox,value_type* arr=NULL)
				
  {
    int nSpareLow[NDIM]={0};
    int nSpareHigh[NDIM]={0};
    gatherSubBoxC(pMin,pMax,nSpareLow,nSpareHigh,localSubBox,arr);
  }

  template <class DT, class IT>
  void gatherSubBoxC(DT pMin[NDIM], DT pMax[NDIM],  
				IT nSpareLow[NDIM], IT nSpareHigh[NDIM], 
				LocalGrid &localSubBox, value_type* arr=NULL)
  {
    IT iMin[NDIM];
    IT iMax[NDIM];
    
    for (int i=0;i<NDIM;++i)
      {
	const std::vector<double> &vCoords = getVertexCoord(i);
	iMin[i]=std::distance(vCoords.begin(),
			      hlp::findFirstHigher(vCoords.begin(),vCoords.end(),pMin[i]));
	iMax[i]=std::distance(vCoords.begin(),
			      hlp::findFirstHigher(vCoords.begin(),vCoords.end(),pMax[i]));
	iMin[i]--;
      }
    
    gatherSubBoxI(iMin,iMax,nSpareLow,nSpareHigh,localSubBox,arr);
  }

  template <class IT>
  void gatherSubBoxI(IT iMin[NDIM], IT iMax[NDIM], 
				LocalGrid &localSubBox,
				value_type* arr=NULL)
  {
    IT nSpareLow[NDIM]={0};
    IT nSpareHigh[NDIM]={0};
    gatherSubBoxC(iMin,iMax,nSpareLow,nSpareHigh,localSubBox,arr);
  }

  // THIS IS NOT WORKING YET !!!
  template <class IT>
  void gatherSubBoxI(IT iMin[NDIM], IT iMax[NDIM], 
				IT nSpareLow[NDIM], IT nSpareHigh[NDIM],
				LocalGrid &localSubBox, value_type* arr=NULL)
  {    
    Params p=params;    
    for (int i=0;i<NDIM;++i)
      {
	IT min = std::max(iMin[i]-nSpareLow[i],static_cast<IT>(0));
	IT max = std::min(iMax[i]+nSpareHigh[i],static_cast<IT>(getResolution(i)));
	p.resolution[i]=max-min;
	p.x0[i]=getVertexCoord(i)[min];
	p.delta[i]=getVertexCoord(i)[max];
	p.delta[i]-=p.x0[i];
	p.lowMargin[i]=p.highMargin[i]=0;   
      }
    
    localSubBox.setName(dataName.c_str());
    localSubBox.initialize(p,arr);
  }
*/
  
  /** \brief Gather from lg the values necessary to apply Kernel at 
   *   coordinates stored in container.
   * \warning Fields must be consecutive or nFields must be 1
   * \warning margins with periodic boundaries (i.e. crossing periodic boundaries) are 
   * not supported !
   */
  template <class Kernel, class Container, bool CellCenteredValues=true>
  void gatherSubsetAtCoords(const Container &container, 
			    LocalGrid &lg, 		    
			    bool alwaysInitializeGrid=false,
			    int nThreads = glb::num_omp_threads,
			    int whichField=-1)
  {
    typedef typename LocalGrid::Data LocalGridData;
    typedef typename Container::iterator Iterator;
    typedef typename Iterator::value_type CoordPtr;
    typedef unsigned int Index;

    //typedef typename hlp::MinimalIntegerType<sizeof(Data)*8,false>::Type IndexType;

    Params localParams = params;
    const long sz=mpiCom->size();
    const long rank=mpiCom->rank();
    int fieldIndex = (whichField<0)?(-whichField-1):whichField;

    // nFields may have been reinterpreted
    localParams.nFields = getNFields();   
    
    if ((alwaysInitializeGrid)||(!lg.isInitialized())) 
      lg.initialize(localParams);
	
    lg.setName(dataName.c_str());    
    
    // Compute the dimensions of the local grid subsets belonging to each MPI process    
    int pos[sz][NDIM];
    int iMax[sz][NDIM];
    int subDim[sz][NDIM];
    long volume[sz];
    for (int rk=0;rk<sz;++rk)
      {
	//receiveCount[i]=1;
	volume[rk]=1;
	for (int j=0;j<NDIM;++j)
	  {
	    //receiveCount[i]*=remoteData[i].params.resolution[j];
	    pos[rk][j]=remoteData[rk].params.position[j]
	      +remoteData[rk].params.lowMargin[j];
	    subDim[rk][j]=remoteData[rk].params.resolution[j]
	      - remoteData[rk].params.lowMargin[j]
	      - remoteData[rk].params.highMargin[j];
	    iMax[rk][j]=pos[rk][j]+subDim[rk][j];
	    volume[rk]*=subDim[rk][j];
	  }
	//receiveDisp[i]=std::distance(getDataPtr(),getDataPtr(pos));
	//printf("RAnk %ld: process %d @[%d %d], size=(%d %d)\n",rank,rk,pos[rk][0],pos[rk][1],subDim[rk][0],subDim[rk][1]);
      }

    // => only 1 MPI process
    if (sz<2)
      {
	if ((getNFields() == 1)||(whichField<0))
	  std::copy_n(getDataPtr(),getLocalNValues()*getNFields(),lg.getDataPtr());
	else
	  {
	    typename LocalGrid::oneFieldIterator itOut(lg.subbox_begin(pos[rank],iMax[rank]),fieldIndex);
	    local_oneFieldIterator itIn(grid->begin(),fieldIndex);
	    const local_oneFieldIterator itIn_end(grid->end(),fieldIndex);
	    std::copy(itIn,itIn_end,itOut);
	  }
	//std::copy_n(getDataPtr(),getLocalNValues()*getNFields(),lg.getDataPtr());
	return;
      }    

    // Reset the local grid
    if (getNFields()==1) std::fill_n(lg.getDataPtr(),lg.getNValues(),0);      
    else
      {
	std::fill(typename LocalGrid::oneFieldIterator(lg.begin(),fieldIndex),
		  typename LocalGrid::oneFieldIterator(lg.end(),fieldIndex),0);
      }

    //std::fill_n(lg.getDataPtr(),lg.getNValues(),0);

    // global grid parameters
    const long * __restrict stride = lg.getValueStride(); 
    long resolution[NDIM];
    for (int i=0;i<NDIM;++i) resolution[i]=lg.getResolution(i);  

    // Set to 1 any pixel that intersects the kernel around any given point coordinate.    
    //gridKernel::BlockStencilT<NDIM,3> kernel;
    Kernel kernel;
    kernel.initialize(resolution,stride,(IS_INTERLEAVED)?getNFields():1);
#pragma omp parallel for num_threads(nThreads) 
    for (int th=0;th<nThreads;++th)
      {	
	long iCoord[NDIM];		

	const auto it_end=container.end(th,nThreads);
	for (auto it=container.begin(th,nThreads); it!=it_end; ++it)
	  {
	    LocalGridData * __restrict ptr=lg.getDataPtr(0,fieldIndex);

	    ptr += kernel.template getCorePixelCoords<IS_PERIODIC,CellCenteredValues>
	      (*it,lg.getOrigin(),lg.getBoxSizeInv(),iCoord);

	    kernel.template apply<IS_PERIODIC,gridKernel::IMPRINT>(ptr,1.0,iCoord);	  

	    /*
	    lg.coordsToICoords(*it,iCoord);	 
	    LocalGridData * __restrict ptr=lg.getDataPtr(0,fieldIndex);	    
	    for (int i=0;i<NDIM;++i) ptr+=iCoord[i]*stride[i];
	    CICkernel.imprint<IS_PERIODIC>(ptr,1.0,iCoord);
	    */

	    /*
	    for (int d=0;d<NDIM;++d)
	      {
		LocalGridData * __restrict p=ptr-stride[d]*bufferSize;
		long count = (bufferSize*2+1)*stride[d];
	
		if (iCoord[d]<bufferSize)
		  {
		    if (IS_PERIODIC)
		      {
			LocalGridData * __restrict p2=p+resolution[d]*stride[d];
			for (int i=0;i<bufferSize-iCoord[d];++i)
			  p2[i*stride[d]] = 1;		
		      }		    
		    p+=(bufferSize-iCoord[d])*stride[d];
		  }

		if (iCoord[d]+bufferSize>=resolution[d]) 
		  {
		    count-=(iCoord[d]+bufferSize-resolution[d]+1)*stride[d];
		    if (IS_PERIODIC)
		      {
			LocalGridData * __restrict p2=p-resolution[d]*stride[d];
			for (int i=count;i<(bufferSize*2+1);++i)
			  p2[i*stride[d]] = 1;	
		      }		    
		  }
		
		for (int i=0;i<count;i+=stride[d]) p[i]=1;
	      }
	    */
	  }
      }

    // Count how many pixels at least must be transfered to each MPI process and compute
    // the dimension of the smallest bounding box containing them all.
    std::vector<int> nSend(sz);
    long bboxVolume[sz];
    int bbox[sz][2][NDIM];
    for (int rk=0;rk<sz;++rk)
      {
	std::copy_n(resolution,NDIM,bbox[rk][0]);
	std::fill_n(bbox[rk][1],NDIM,0);
      }
	
#pragma omp parallel for num_threads(nThreads)
    for (int rk=0;rk<sz;++rk)
      {
	if (rk==mpiCom->rank()) 
	  {
	    bboxVolume[rk]=1;
	    for (int i=0;i<NDIM;++i)
	      {
		bbox[rk][0][i]=0;
		bbox[rk][1][i]=subDim[rk][i]-1;
		bboxVolume[rk]*=subDim[rk][i];
	      }
	    nSend[rk]=bboxVolume[rk];
	    continue;
	  }

	typename LocalGrid::oneFieldIterator it(lg.subbox_begin(pos[rk],iMax[rk]),fieldIndex);
	const typename LocalGrid::oneFieldIterator it_end(lg.subbox_end(pos[rk],iMax[rk]),fieldIndex);
	nSend[rk]=0;
	for (;it!=it_end;++it)
	  {
	    if (*it == 1)
	      {
		nSend[rk]++;
		const int *w=it.get_w();
		for (int i=0;i<NDIM;++i)
		  {
		    if (w[i]<bbox[rk][0][i]) bbox[rk][0][i]=w[i];
		    if (w[i]>bbox[rk][1][i]) bbox[rk][1][i]=w[i];
		  }
	      }
	  }

	bboxVolume[rk]=1;
	for (int i=0;i<NDIM;++i) 
	  bboxVolume[rk]*=(bbox[rk][1][i]-bbox[rk][0][i]+1);
      }

    // Correct bboxVolume when the intersection is void   
#pragma omp parallel for num_threads(nThreads)
    for (int rk=0;rk<sz;++rk)
      {
	if (nSend[rk]==0)
	  bboxVolume[rk]=0;	
      }

    /*
    long nSendTotal[3]={0};
    for (int rk=0;rk<sz;++rk)
      {
	if (rk!=mpiCom->rank())
	  {
	    nSendTotal[0]+=nSend[rk];
	    nSendTotal[1]+=bboxVolume[rk];
	    nSendTotal[2]+=volume[rk];
	  }
      }
    glb::console->print<LOG_STD_ALL>("-> Exchanging items: %ld indexed / %ld boxed / %ld all.\n",
				     nSendTotal[0],nSendTotal[1],nSendTotal[2]);    
    */

    // Select the exchange mode and update nSend accordingly
    int exchangeMode[sz][2]; //[i][0]=>send mode, [i][1]=>receive mode
    int nSendCum[sz+1]; 
    
    //double ratio = double(sizeof(Data))/double(sizeof(Index));
    // Both modes exchange the same amount of info if a 'threshold' 
    // fraction of the voxels has to be retrieved
    double threshold=double(sizeof(Data))/double(sizeof(Index)+sizeof(Data)); //ratio/(1.0+ratio); 
    nSendCum[0]=0;
    for (int rk=0;rk<sz;++rk)
      {
	// store nSend in exchangeMode[][0] temporarily for MPI convenience
	if (nSend[rk]>(volume[rk]*threshold))
	  {
	    nSend[rk]=0;
	    
	    exchangeMode[rk][0]=nSend[rk]; 
	    exchangeMode[rk][1]=1; // send nothing, receive everything
	  }
	else 
	  {	    
	    exchangeMode[rk][0]=nSend[rk];
	    exchangeMode[rk][1]=0; // send indexes, receive subset
	  }

	nSendCum[rk+1]=nSendCum[rk]+nSend[rk];
      }
 
    std::vector<int> nReceive(sz);
    int nReceiveCum[sz+1];
    int buffer[sz][2];

    if (glb::console->willPrint<LOG_DEBUG>())
      {
	std::stringstream ssR;
	
	ssR << "Rank "<<rank<<" : ";
	for (int rk=0;rk<sz;++rk) 
	  ssR << "("<<rk<<","<<exchangeMode[rk][1]<<":"<<exchangeMode[rk][0]<<")";

	glb::console->print<LOG_DEBUG>("%s\n",ssR.str().c_str());
      }

    // Exchange the mode and how many pixels to send
    mpiCom->Alltoall(&exchangeMode[0][0],&buffer[0][0],2);

    for (int rk=0;rk<sz;++rk)
      {
	exchangeMode[rk][0]=buffer[rk][1]; // receive mode
	nReceive[rk]=buffer[rk][0];
      }

    //mpiCom->Alltoall(nSend,nReceive);

    nReceiveCum[0]=0;
    for (int rk=0;rk<sz;++rk)
      nReceiveCum[rk+1]=nReceiveCum[rk]+nReceive[rk];
    
    /*
    size_t nSendAlloc = std::max(sizeof(Index)*nSendCum[sz],
				 sizeof(Data)*nSendDataCum[sz]);
    size_t nReceiveAlloc = std::max(sizeof(Index)*nReceiveCum[sz],
				    sizeof(Data)*nReceiveDataCum[sz]);
    */
    Index *sendBuffer = (Index*)malloc(sizeof(Index)*nSendCum[sz]);
    Index *receiveBuffer = (Index*)malloc(sizeof(Index)*nReceiveCum[sz]);
    
    //std::vector<Index> sendBuffer(nSendCum[sz]);
    //Compute the indices to exchange if necessary
#pragma omp parallel for num_threads(nThreads)
    for (int rk=0;rk<sz;++rk)
      {
	if (nSend[rk]==0) continue;

	typename LocalGrid::oneFieldIterator it(lg.subbox_begin(pos[rk],iMax[rk]),fieldIndex);
	const typename LocalGrid::oneFieldIterator it_end(lg.subbox_end(pos[rk],iMax[rk]),fieldIndex);
	Index index=0;
	Index *indexVec = &sendBuffer[nSendCum[rk]];

	for (;it!=it_end;++it,++index)
	  {
	    if (*it == 1)
	      {
		(*indexVec)=index;
		indexVec++;
	      }
	  }	
      }

    if (glb::console->willPrint<LOG_DEBUG>())
      {
	std::stringstream ssS;
	std::stringstream ssR;
	
	ssS << "Rank "<<rank<<" sends : ";
	for (int rk=0;rk<sz;++rk) 
	  if (nSend[rk]) ssS << "("<<rk<<","<<nSend[rk]<<")";
	
	ssR << "Rank "<<rank<<" recvs : ";
	for (int rk=0;rk<sz;++rk) 
	  if (nReceive[rk]) ssR << "("<<rk<<","<<nReceive[rk]<<")";
	glb::console->print<LOG_DEBUG>("%s\n",ssS.str().c_str());
	glb::console->print<LOG_DEBUG>("%s\n",ssR.str().c_str());
      }
    
    
    // exchange the indices of requested pixels
    //std::vector<Index> receiveBuffer(nReceiveCum[sz]);
    mpiCom->Alltoallv(&sendBuffer[0],&nSend[0],&nSendCum[0],
		      &receiveBuffer[0],&nReceive[0],&nReceiveCum[0]);
 
    // Just creating aliases for readability ...
    int *nSendData=&nReceive[0];
    int *nSendDataCum=&nReceiveCum[0];
    int *nReceiveData=&nSend[0];
    int *nReceiveDataCum=&nSendCum[0];
    
    // allocate buffers to store pixel values
    free(sendBuffer);sendBuffer=NULL;
    Data *sendDataBuffer = (Data*) malloc(sizeof(Data)*nSendDataCum[sz]);    
    
    // copy the value of requested pixel
    // FIXME: margins with PBC wont work here !
#pragma omp parallel for num_threads(nThreads)
    for (int rk=0;rk<sz;++rk)
      {
	if (nSendData[rk]==0) continue;

	Data *dataVec = &sendDataBuffer[nSendDataCum[rk]];
	Index *indexVec = receiveBuffer + nReceiveCum[rk];
	Data *ptr=getDataPtr(0,fieldIndex);

	if (IS_INTERLEAVED)
	  {
	    long fac = getNFields();
	    for (Index i=0;i<nReceive[rk];++i)
	      dataVec[i] = ptr[fac*indexVec[i]];
	  }
	else
	  {
	    for (Index i=0;i<nReceive[rk];++i)
	      dataVec[i] = ptr[indexVec[i]];
	  }	
      }
      
    if (glb::console->willPrint<LOG_DEBUG>())
      {
	std::stringstream ssS;
	std::stringstream ssR;
	ssS << "Rank "<<rank<<" sends: ";
	for (int rk=0;rk<sz;++rk) 
	  if (nSendData[rk]) ssS << "("<<rk<<","<<nSendData[rk]<<")";
	ssR << "Rank "<<rank<<" recvs: ";
	for (int rk=0;rk<sz;++rk) 
	  if (nReceiveData[rk]) ssR << "("<<rk<<","<<nReceiveData[rk]<<")";
	glb::console->print<LOG_DEBUG>("%s\n",ssS.str().c_str());
	glb::console->print<LOG_DEBUG>("%s\n",ssR.str().c_str());
      }
    
    // and finally exchange requested pixels
    free(receiveBuffer);receiveBuffer=NULL;
    Data *receiveDataBuffer = (Data*) malloc(sizeof(Data)*nReceiveDataCum[sz]);
    mpiCom->Alltoallv(&sendDataBuffer[0],&nSendData[0],&nSendDataCum[0],
		      &receiveDataBuffer[0],&nReceiveData[0],&nReceiveDataCum[0]);
    
#pragma omp parallel for num_threads(nThreads)
    for (int rk=0;rk<sz;++rk)
      {
	if (nReceiveData[rk]==0) continue;

	typename LocalGrid::oneFieldIterator it(lg.subbox_begin(pos[rk],iMax[rk]),fieldIndex);
	Data *dataVec = &receiveDataBuffer[nReceiveDataCum[rk]];
	int nFound=0;

	do {
	  if (*it == 1) (*it)=dataVec[nFound++];
	  ++it;
	} while(nFound<nReceiveData[rk]);
      }

    free(sendDataBuffer);
    free(receiveDataBuffer); 
    
    // Now we still have to exchange mode 1 pixels 
    // (i.e sending everything)

    /*
    int old_nSendDataCum=nSendDataCum[sz];
    int old_nReceiveDataCum=nReceiveDataCum[sz];
    nSendCum[0]=0;nReceiveCum[0]=0;
    for (int rk=0;rk<sz;++rk)
      {
	if (rk==mpiCom->rank())
	  {
	    nReceiveData[rk]=0;
	    nSendData[rk]=0;
	  }

	if (exchangeMode[rk][1]==1)
	  nReceiveData[rk]=volume[rk];
	else 
	  nReceiveData[rk]=0;

	if (exchangeMode[rk][1]==1)
	  nSendData[rk]=volume[rk];
	else 
	  nSendData[rk]=0;

	nSendDataCum[rk+1]=nSendDataCum[rk]+nSendData[rk];
	nReceiveDataCum[rk+1]=nReceiveDataCum[rk]+nReceiveData[rk];
      }
    
    if (nReceiveDataCum[sz]>old_nReceiveDataCum)
      {
	free(receiveDataBuffer);receiveDataBuffer=NULL;
	Data *receiveDataBuffer = (Data*) malloc(sizeof(Data)*nReceiveDataCum[sz]);
      }
    */
    /*
    if (nSendDataCum[sz]>old_nSendDataCum)
      {
	free(sendDataBuffer);sendDataBuffer=NULL;
	Data *sendDataBuffer = (Data*) malloc(sizeof(Data)*nSendDataCum[sz]);
      }
    */
    /*
    mpiCom->Alltoallv(&sendDataBuffer[0],&nSendData[0],&nSendDataCum[0],
		      &receiveDataBuffer[0],&nReceiveData[0],&nReceiveDataCum[0]);
    
    if (nReceiveDataCum[sz]>0)
      {
	for (int rk=0;rk<sz;++rk)
	  {
	    if (nReceiveData[rk])
	      {
		typename LocalGrid::oneFieldIterator itOut(lg.subbox_begin(pos[rk],iMax[rk]),fieldIndex); 
		std::copy_n(receiveDataBuffer[rk],nReceiveData[rk],itOut);
	      }
	  }
      }
    */    

    // Finally copy the locally stored part of the grid
    /*
    typename LocalGrid::oneFieldIterator itOut(lg.subbox_begin(pos[rank],iMax[rank]),fieldIndex); 
    local_oneFieldIterator itIn(grid->begin(),fieldIndex);
    const local_oneFieldIterator itIn_end(grid->end(),fieldIndex);
    std::copy(itIn,itIn_end,itOut);   
    */

    //free(sendDataBuffer);
    
    // Gather the next field if necessary (not optimal but it works ...)
    if ((whichField<0)&&(fieldIndex<getNFields()-1))
      gatherSubsetAtCoords<Kernel>(container,lg,false,nThreads,whichField-1);
    else 
      {
	// We still have to transfer mode 1 chunks and local piece
	gatherAll(lg,exchangeMode,false);
      }
  }

  // sendReceive[rank][0] is true if sending to process 'rank'
  // sendReceive[rank][1] is true if receiving from process 'rank'
  // Note that it is the user's responsibility to give consistant values, 
  // otherwise your should rather use the slower version where only 
  // receiveFrom is specified
  template <class IT>
  void gatherAll(LocalGrid &lg,const IT sendReceive[][2], bool alwaysInitializeGrid=false)
  {
    Params localParams = params;
    const long sz=mpiCom->size();

    // because nfields may have been reinterpreted !!!
    localParams.nFields = getNFields();   
    
    if ((alwaysInitializeGrid)||(!lg.isInitialized())) 
      lg.initialize(localParams);
	
    lg.setName(dataName.c_str());

    if (sz<2)
      {
	std::copy_n(getDataPtr(),getLocalNValues()*getNFields(),lg.getDataPtr());
	return;
      }

    for (int i=0;i<NDIM;++i)
      {
	if ((grid->getLowMargin(i)!=0)||(grid->getHighMargin(i)!=0))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("gatherAll is not compatible with margins yet.\n");
	    glb::console->print<LOG_ERROR>("=> this can be fixed by adapting the sendTypes !\n");
	    exit(-1);
	  }
      }
    
    int pos[sz][NDIM];
    int subDim[sz][NDIM];    
    for (int i=0;i<sz;++i)
      {
	for (int j=0;j<NDIM;++j)
	  {
	    pos[i][j]=remoteData[i].params.position[j] + remoteData[i].params.lowMargin[j];
	    subDim[i][j]=remoteData[i].params.resolution[j] - 
	      remoteData[i].params.lowMargin[j] -
	      remoteData[i].params.highMargin[j];
	  }
      }

    int res[NDIM];    
    std::copy_n(params.resolution,NDIM,res);

    MPI_Datatype sendType[sz];
    MPI_Datatype receiveType[sz];

    long nSend = getLocalNValues();
    std::vector<int> sendCount(sz,nSend);//[sz]={nSend};
    std::vector<int> sendDisp(sz,0);
    std::vector<int> receiveDisp(sz,0);
    std::vector<int> receiveCount(sz,1);

    for (int rk=0;rk<sz;++rk)
      {
	if (!sendReceive[rk][0]) sendCount[rk]=0;
	if (!sendReceive[rk][1]) receiveCount[rk]=0;
      }
   
    if (getNFields()<2)
      {	
	for (int i=0;i<sz;++i)
	  {
	    sendType[i]=MPI_Type<Data>::get();
	    if (receiveCount[i]==0) receiveType[i]=MPI_Type<Data>::get();
	    else {
	    MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
				     MPI_Type<Data>::get(),&receiveType[i]);
	    MPI_Type_commit(&receiveType[i]);
	    }
	  }

	mpiCom->Alltoallw(getDataPtr(),&sendCount[0],&sendDisp[0],sendType,
			  lg.getDataPtr(),&receiveCount[0],&receiveDisp[0],receiveType);
      }
    else
      {
	if (IS_INTERLEAVED)
	  {
	    MPI_Datatype fieldType;
	    MPI_Type_contiguous(getNFields(),MPI_Type<Data>::get(),&fieldType);
	    MPI_Type_commit(&fieldType);

	    for (int i=0;i<sz;++i)
	      {
		sendType[i]=fieldType;
		if (receiveCount[i]==0) receiveType[i]=MPI_Type<Data>::get();
		else {
		MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
					 fieldType,&receiveType[i]);
		MPI_Type_commit(&receiveType[i]);
		}	
	      }

	    mpiCom->Alltoallw(getDataPtr(),&sendCount[0],&sendDisp[0],sendType,
			      lg.getDataPtr(),&receiveCount[0],&receiveDisp[0],receiveType);

	    MPI_Type_free(&fieldType);
	  }
	else
	  {
	    for (int i=0;i<sz;++i)
	      {
		sendType[i]=MPI_Type<Data>::get();
		if (receiveCount[i]==0) receiveType[i]=MPI_Type<Data>::get();
		else {
		MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
					 MPI_Type<Data>::get(),&receiveType[i]);
		MPI_Type_commit(&receiveType[i]);
		}
	      }

	    for (int i=0;i<getNFields();++i)
	      {
		mpiCom->Alltoallw(getDataPtr(0,i),
				  &sendCount[0],&sendDisp[0],sendType,
				  lg.getDataPtr(0,i),
				  &receiveCount[0],&receiveDisp[0],receiveType);
	      }
	  }
      }

    for (int i=0;i<sz;++i) 
      if (receiveCount[i]!=0) 
	MPI_Type_free(&receiveType[i]);
  }

  void gatherAll(LocalGrid &lg, bool alwaysInitializeGrid=false)
  {
    int size = mpiCom->size();
    char sendReceive[size][2];
    std::fill_n(&sendReceive[0][0],size*2,1);

    gatherAll(lg,sendReceive,alwaysInitializeGrid);
  }
  
  template <class T>
  void gatherAll(LocalGrid &lg,const T receiveFrom[], bool alwaysInitializeGrid=false)
  {
    int size = mpiCom->size();
    T sendTo[size];
    T sendReceive[size][2];
    
    mpiCom->Alltoall(receiveFrom,sendTo,1);
    for (int rk=0;rk<size;++rk)
      {
	sendReceive[rk][0]=sendTo[rk];
	sendReceive[rk][1]=receiveFrom[rk];
      }

    gatherAll(lg,sendReceive,alwaysInitializeGrid);
  }

  /*
  //void gatherAll(LocalGrid &lg, int target=-1, bool alwaysInitializeGrid=false)
  // This does not work when there are margins (have to adapt the sendtypes !)
  void gatherAll(LocalGrid &lg, bool alwaysInitializeGrid=false)
  {
    Params localParams = params;
    //bool isTarget = ((target<0)||(target==mpiCom->rank()));
    long sz=mpiCom->size();

    // because nfields may have been reinterpreted !!!
    localParams.nFields = getNFields();   
    
    if ((alwaysInitializeGrid)||(!lg.isInitialized())) 
      lg.initialize(localParams);
	
    lg.setName(dataName.c_str());

    if (sz<2)
      {
	std::copy_n(getDataPtr(),getLocalNValues()*getNFields(),lg.getDataPtr());
	return;
      }

    for (int i=0;i<NDIM;++i)
      {
	if ((grid->getLowMargin(i)!=0)||(grid->getHighMargin(i)!=0))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("gatherAll is not compatible with margins yet.\n");
	    glb::console->print<LOG_ERROR>("=> this can be fixed by adapting the sendTypes !\n");
	    exit(-1);
	  }
      }
    
    //int receiveCount[sz];
    //int receiveDisp[sz];
    int pos[sz][NDIM];
    int subDim[sz][NDIM];    
    for (int i=0;i<sz;++i)
      {
	//receiveCount[i]=1;
	for (int j=0;j<NDIM;++j)
	  {
	    //receiveCount[i]*=remoteData[i].params.resolution[j];
	    pos[i][j]=remoteData[i].params.position[j] + remoteData[i].params.lowMargin[j];
	    subDim[i][j]=remoteData[i].params.resolution[j] - 
	      remoteData[i].params.lowMargin[j] -
	      remoteData[i].params.highMargin[j];
	  }
	//receiveDisp[i]=std::distance(getDataPtr(),getDataPtr(pos));
      }

    int res[NDIM];    
    std::copy_n(params.resolution,NDIM,res);

    MPI_Datatype sendType[sz];
    MPI_Datatype receiveType[sz];

    long nSend = getLocalNValues();
    std::vector<int> sendCount(sz,nSend);//[sz]={nSend};
    std::vector<int> sendDisp(sz,0);
    std::vector<int> receiveDisp(sz,0);
    std::vector<int> receiveCount(sz,1);
    // int sendCount[sz]={nSend};
    // int sendDisp[sz]={0};
    // int receiveDisp[sz]={0};
    // int receiveCount[sz]={1};

    if (getNFields()<2)
      {	
	for (int i=0;i<sz;++i)
	  {
	    sendType[i]=MPI_Type<Data>::get();
	    MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
				     MPI_Type<Data>::get(),&receiveType[i]);
	    MPI_Type_commit(&receiveType[i]);
	  }

	mpiCom->Alltoallw(getDataPtr(),&sendCount[0],&sendDisp[0],sendType,
			  lg.getDataPtr(),&receiveCount[0],&receiveDisp[0],receiveType);
      }
    else
      {
	if (IS_INTERLEAVED)
	  {
	    MPI_Datatype fieldType;
	    MPI_Type_contiguous(getNFields(),MPI_Type<Data>::get(),&fieldType);
	    MPI_Type_commit(&fieldType);

	    for (int i=0;i<sz;++i)
	      {
		sendType[i]=fieldType;
		MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
					 fieldType,&receiveType[i]);
		MPI_Type_commit(&receiveType[i]);
	      }

	    mpiCom->Alltoallw(getDataPtr(),&sendCount[0],&sendDisp[0],sendType,
			      lg.getDataPtr(),&receiveCount[0],&receiveDisp[0],receiveType);

	    MPI_Type_free(&fieldType);
	  }
	else
	  {
	    for (int i=0;i<sz;++i)
	      {
		sendType[i]=MPI_Type<Data>::get();
		MPI_Type_create_subarray(NDIM,res,subDim[i],pos[i],MPI_ORDER_FORTRAN,
					 MPI_Type<Data>::get(),&receiveType[i]);
		MPI_Type_commit(&receiveType[i]);
	      }

	    for (int i=0;i<getNFields();++i)
	      {
		mpiCom->Alltoallw(getDataPtr(0,i),
				  &sendCount[0],&sendDisp[0],sendType,
				  lg.getDataPtr(0,i),
				  &receiveCount[0],&receiveDisp[0],receiveType);
	      }
	  }
      }

    for (int i=0;i<sz;++i) 
      MPI_Type_free(&receiveType[i]);
  }
  */
  
  template <class F>
  double quadrature(F &f, int field=0,
		    int nThreads=glb::num_omp_threads)
  {
    double result=getLocalGrid()->quadrature(f,field,nThreads);
    return mpiCom->sum(result);
  }

private:

  template <class AMR, class T>
  void addAmrGrid_getParams(int index, T pos[NDIM], T resolution[NDIM], T &level)
  {
    const Params &rParams = remoteData[index].params;
    level = AMR::MAX_LEVEL;

    long levelResolution = 1L<<level;
    long parentResolution[NDIM];
    

    if (rParams.haveParentGrid)
      {
	std::copy(rParams.position,rParams.position+NDIM,pos);
	std::copy(rParams.parentResolution,rParams.parentResolution+NDIM,parentResolution);
      }
    else
      {
	std::fill(pos,pos+NDIM,0);
	std::copy(rParams.resolution,rParams.resolution+NDIM,parentResolution);
      }

    std::copy(rParams.resolution, rParams.resolution+NDIM,resolution);
      
    for (int i=0;i<NDIM;++i)
      {
	long factor = levelResolution/parentResolution[i];
	if (factor*parentResolution[i] != levelResolution)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("This function is only implemented for grids with sizes a power of 2.\n");
	    exit(-1);
	  }
	pos[i]*=factor;
	resolution[i]*=factor;
      }
    level -= AMR::ROOT_LEVEL;
  }

  struct RemoteData
  {
    Params params;
    //int gridPos[NDIM];

    template <class W>
    void write(W *writer)
    {
      params.write(writer);
      //writer->write(gridPos,NDIM);
    }

    template <class R>
    void read(R *reader)
    {
      params.read(reader);
      //reader->read(gridPos,NDIM);
    }

    template <class Slicer>
    void set(const Slicer &slicer)
    {
      params=slicer.getLocalGridParams();
      //for (int j=0;j<NDIM;++j)
      //gridPos[j]=slicer.getLocalGridPos(j);
    }
  };

  Params params;
  int numThreads;
  
  MpiCommunication *mpiCom;
  long nGrids;
  LocalGrid *grid;

  std::string dataName;

  double xMin[NDIM];
  double xMax[NDIM];
  int sliceCount[NDIM];
  int slicePos[NDIM];
  //int localGridPos[NDIM];

  //int neighborIndex[NDIM][2];

  std::vector<double> vertexCoord[NDIM];
  std::vector<double> cellCoord[NDIM];
  std::vector<double> valueCoord[NDIM];

  int onLowBoundary[NDIM];
  int onHighBoundary[NDIM];
  int onBoundary[NDIM];

  int initialized;

  //std::vector<Params> allParams;
  //std::vector<int> remoteGridPos[NDIM];
  std::vector<RemoteData> remoteData;
  GridTopology topology;
};

/** \}*/
#include "../internal/namespace.footer"
#endif
