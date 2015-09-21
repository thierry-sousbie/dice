#ifndef __REGULAR_GRID_SLICER_SIMPLE_HXX__
#define __REGULAR_GRID_SLICER_SIMPLE_HXX__

#include "../dice_globals.hxx"

#include "./internal/regularGridSlicerBase.hxx"

#include "./valLocationType.hxx"

#include "../tools/MPI/mpiCommunication.hxx"

/**
 * @file 
 * @brief A simple MPI slicer for a global regular grid 
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup GRID
 *   \{
 */

/**
 * \class RegularGridSlicerSimpleT
 * \brief A simple grid slicer that cuts a global grid into slabs of roughly equal size
 * that may overlap
 * \tparam G A regular grid class
 */
template <class G>
class RegularGridSlicerSimpleT : 
  public internal::RegularGridSlicerBaseT<G> 
{
public:
  typedef internal::RegularGridSlicerBaseT<G> Base;
  typedef RegularGridSlicerSimpleT<G> MyType;

  static const int NDIM = G::NDIM;

  typedef G Grid;
  typedef typename Grid::Params Params;
  typedef typename Grid::Scale Scale;
  //typedef typename Grid::ValLocationType ValLocationType;
  //typedef typename Grid::ValLocationTypeV ValLocationTypeV;
  typedef typename Scale::ScaleTypeV ScaleTypeV;
  typedef typename Scale::ScaleType  ScaleType;
  typedef typename Grid::GridNav GridNav;
  typedef typename Grid::Direction Direction;

  template <class I>
  RegularGridSlicerSimpleT(const Params &gp_, long nChunks_,
			   MpiCommunication *com_, 
			   int nThreads_,
			   const I lowMarginSize_[NDIM],
			   const I highMarginSize_[NDIM],
			   bool periodic_=false)
  {
    initialize(gp_,nChunks_,com_,nThreads_,lowMarginSize_,highMarginSize_,periodic_);   
  }

  RegularGridSlicerSimpleT(const Params &gp_, long nChunks_,
			   MpiCommunication *com_, 
			   int nThreads_,
			   bool periodic_=false)
  {
    long lowMarginSize_[NDIM]={0};
    long highMarginSize_[NDIM]={0};
    initialize(gp_,nChunks_,com_,nThreads_,lowMarginSize_,highMarginSize_,periodic_);
  }
  /*
  RegularGridSlicerSimpleT()
  {}
  */
  void slice(int index)
  {    
    neighbors.clear();    
    long nSlices = nChunks;
    globalGridParams=gp;
    sliceCount[0]=nSlices;
    for (int i=1;i<NDIM;++i) sliceCount[i]=1;
    myIndex=index;
    slicePos[0]=index;
    for (int i=1;i<NDIM;++i) slicePos[i]=0;   
    
    int n0,n1;
    this->getChunkIndices(gp,slicePos[0],sliceCount[0],0,n0,n1);
    gridParams=this->divide(gp,n0,n1,lowMarginSize[0],highMarginSize[0],periodic,0);
			    
    /*
    gridParams=this->divide(gp,slicePos[0],sliceCount[0],
			    lowMarginSize[0],highMarginSize[0],
			    periodic,0);
    */
    if (slicePos[0]+1 >= sliceCount[0]) 
      {if (periodic) neighbors.push_back(Neighbor(0,1,0));}
    else neighbors.push_back(Neighbor(0,1,myIndex+1));
    
    if (slicePos[0]<=0)
      {if (periodic) neighbors.push_back(Neighbor(0,-1,sliceCount[0]-1));}
    else neighbors.push_back(Neighbor(0,-1,myIndex-1));

    gridParams.haveParentGrid=true;
    gridParams.parentNDim = NDIM;    
    gridParams.minElementsCount=0;    

    for (int i=1;i<NDIM;++i) gridParams.position[i]=0;
    for (int i=0;i<NDIM;++i)
      {
	gridParams.parentResolution[i]= gp.resolution[i];
	gridParams.parentX0[i] = gp.x0[i];
	gridParams.parentDelta[i] = gp.delta[i];
      }    
  }

  const Params &getLocalGridParams() const
  {
    return gridParams;
  }

  const Params &getGlobalGridParams() const
  {
    return globalGridParams;
  }
  
  long neighborsCount() const
  {
    return neighbors.size();
  }

  template <class T, class T2>
  void getNeighborInfo(int which, T &neiDim, T &neiDir, T2 &neiIndex) const
  {
    neiDim=neighbors[which].dim;
    neiDir=neighbors[which].dir;
    neiIndex=neighbors[which].index;
  }

  long getSliceCount(int dim) const
  {
    return sliceCount[dim];
  }

  long getSlicePos(int dim) const
  {
    return slicePos[dim];
  }

  MpiCommunication *getMpiCom() const
  {
    return mpiCom;
  }

  long getNChunks()
  {
    return nChunks;
  }

  long getNThreads()
  {
    return nChunks;
  }

private:
  template <typename I>
  void initialize(const Params &gp_, long nChunks_,
		  MpiCommunication *com_, 
		  int nThreads_,
		  const I lowMarginSize_[NDIM],
		  const I highMarginSize_[NDIM],
		  bool periodic_)
  {
    gp=gp_;
    nChunks=nChunks_;
    mpiCom=com_;
    numThreads=nThreads_;
    std::copy(lowMarginSize_,lowMarginSize_+NDIM,lowMarginSize);
    std::copy(highMarginSize_,highMarginSize_+NDIM,highMarginSize);
    periodic = periodic_; 

     if (nChunks > gp.resolution[0])
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Cannot produce more slices than the size of the first dimension: max is %d, %ld required.\n",gp.resolution[0],nChunks);
	exit(-1);
      }

     slice(mpiCom->rank());
  }

  struct Neighbor
  {
    Neighbor(int dm, int dr, int id):
      dim(dm),dir(dr),index(id)
    {}
   
    void set(int dm, int dr, int id)
    {dim=dm;dir=dr;index=id;}

    int dim;
    int dir;
    int index;
  };

  Params globalGridParams;
  Params gridParams; 

  std::vector< Neighbor > neighbors;
      
  int sliceCount[NDIM];//int slicerSize[DIMS];  // number of slices in each dim
  int myIndex;//int myID;     
  int slicePos[NDIM];//int myCoords[DIMS]; // coords within the slices   
  bool periodic;

  long lowMarginSize[NDIM];
  long highMarginSize[NDIM];
  int numThreads;
  MpiCommunication *mpiCom;
  Params gp;
  long nChunks;
};

/** \}*/
#include "../internal/namespace.footer"
#endif
