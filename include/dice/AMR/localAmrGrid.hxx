#ifndef __LOCAL_AMR_GRID_HXX__
#define __LOCAL_AMR_GRID_HXX__

#include <string>  
#include <iostream> 
#include <sstream>  


#ifdef HAVE_BOOST
#include "../tools/wrappers/boostMultiprecisionFloat128.hxx"
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif

#include "../dice_globals.hxx"

#include "../tools/memory/memoryPool.hxx"
#include "../tools/memory/iterableMemoryPool.hxx"
#include "../tools/helpers/helpers.hxx"

#include "../geometry/boundaryType.hxx"
#include "../geometry/geometricProperties.hxx"

#include "../AMR/localAmrGridVoxel.hxx"
#include "../AMR/localAmrGridRaytracer.hxx"
#include "../AMR/localAmrGridProjector.hxx"
#include "../AMR/localAmrGridVisitors.hxx"

#include "../IO/vtkAmrWriter.hxx"

#include "./internal/localAmrGridDefaultVisitors.hxx"
#include "./internal/localAmrGridDefaultOverlapVisitors.hxx"

/**
 * @file 
 * @brief  A local AMR grid designed to allow the projection of an unstructured mesh
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup AMR
 *   \{
 */

/**
 * \class LocalAmrGridT
 * \brief  A local AMR grid designed to allow the projection of an unstructured mesh
 *
 * Except for template parameters, levels are measured from the root nodes (i.e. a root
 * voxel has level 0 and its first refinement level 1)
 * ML is presently limited to 31 in 2D and 20 in 3D, but this could be easily lifted at 
 * the expanse of some memory ...
 *
 * Point location in vertices for degenerate cases uses simulation of simplicity such that
 * if Xi is any vertex corner coordinate along dimension i (0<=i<NDIM), then
 *            Xi -> Xi - pow(eps,NDIM-i)
 *
 * \tparam ND Number of dimensions
 * \tparam DT Type of data stored within each voxel
 * \tparam BT Boundary conditions type
 * \tparam RL Root level (the minimal level of a voxel)
 * \tparam ML Maximum level (the maximum level a voxel can reach)
 * \tparam MAX_TH Maximum number of threads (default=64) 
 *
 * \todo improve neighbor finding : check N level above the current level if the neighbor
 * is not in a parent cell and start the search from that level if it is the case (if N=4,
 * it will be the case in 90%+ of the cases).
 * \todo implement a hash map version ?
 */

template <long ND, 
	  typename DT=double, 
	  long BT = BoundaryType::NONE, 	  
	  long RL = 6, 
	  long ML = (63/ND-1),
	  int MAX_TH = 64>
class LocalAmrGridT 
{
public:
  typedef LocalAmrGridT<ND,DT,BT,RL,ML> MyType;
  
  typedef DT Data;
  //! ICoord needs to be at least an (ML+1)*ND bytes large unsigned integer type
  typedef typename hlp::MinimalIntegerType<(ML+1)*ND,false>::Type ICoord; 
  //typedef unsigned long ICoord;

  typedef LocalAmrGridVoxelT<MyType> Voxel;
  typedef GeometricPropertiesT<double,ND,ND,BT,BoundaryType::NONE> GeometricProperties;

  typedef LocalAmrGridRaytracerT<MyType> Raytracer;

  template <class U>
  friend class LocalAmrGridVoxelT;
  
  // template <int,class,class,template <class> class> 
  // friend class internal::LocalAmrGridProjectorT;
  // template <int,class,class,template <class> class> 
  // friend class internal::LocalAmrGridProjectorBaseT;
  
  static const long BOUNDARY_TYPE = BT;
  static const long PERIODIC_BOUNDARIES = (BT==BoundaryType::PERIODIC);
  static const long NDIM = ND;   
  //static const double NDIM_INV = 1.0/NDIM;
  static const long ROOT_LEVEL = RL; 
  static const long MAX_LEVEL = ML;
  static const long MAX_LEVEL_FROM_ROOT = ML-RL;
  static const long MAX_LEVEL_FROM_ROOT_PLUS_ONE = ML-RL+1L;
  static const ICoord BBOX_ILEN = (1L<<(MAX_LEVEL+1L));
  static const ICoord ROOT_VOXEL_ILEN = (1L<<(MAX_LEVEL_FROM_ROOT+1L));
  static const ICoord ROOT_VOXEL_IHALFLEN = (1L<<(MAX_LEVEL_FROM_ROOT));
  static const ICoord N_ROOT_PER_DIM = (1L<<ROOT_LEVEL);

  static const long  CHILDREN_COUNT = (1L<<ND);
  static const long ROOT_VOXELS_COUNT = hlp::IntPower<(1L<<ROOT_LEVEL),ND>::Result::value;

  static const ICoord INDEX_DEC  = (MAX_LEVEL+1L);        
  static const ICoord INDEX_MASK = (1L<<(MAX_LEVEL+1L))-1L; // == (1<<INDEX_DEC)-1
  
  /** \brief Define a direction and orientation in space. Directions can be combined using 
   *  the '|' operator.
   * \warning Not all combinations make sense, for instance \a (Left|Rear|Top) is valid 
   * but \a (Left|Right) will be interpreted as \a Right ...
   */
  enum Direction {
    Left   = (0<<16) | (1<<0), /**< Negative orientation along axis 0*/
    Right  = (1<<16) | (1<<0), /**< Positive orientation along axis 0*/

    Rear   = (0<<17) | (1<<1), /**< Negative orientation along axis 1*/
    Front  = (1<<17) | (1<<1), /**< Positive orientation along axis 1*/

    Bottom = (0<<18) | (1<<2), /**< Negative orientation along axis 2*/
    Top    = (1<<18) | (1<<2), /**< Positive orientation along axis 2*/

    Left3  = (0<<19) | (1<<3), /**< Negative orientation along axis 3*/
    Right3 = (1<<19) | (1<<3), /**< Positive orientation along axis 3*/

    Left4  = (0<<20) | (1<<4), /**< Negative orientation along axis 4*/
    Right4 = (1<<20) | (1<<4), /**< Positive orientation along axis 4*/

    Left5  = (0<<21) | (1<<5), /**< Negative orientation along axis 5*/
    Right5 = (1<<21) | (1<<5), /**< Positive orientation along axis 5*/
  };  

  static std::string classHeader() {return "local_amr_grid";}
  static float classVersion() {return 0.10;}
  /** \brief Constructor with bounding box definition
   *  \param x0 an iterator to the lower left coordinates of the bounding box
   *  \param delta an iterator to the size of the bounding box along each dimension
   */
  template <class InputIterator>
  LocalAmrGridT(InputIterator x0, InputIterator delta): 
		//,MpiCommunication *mpiCom_=glb::mpiComWorld):
    //rootVoxels(ROOT_VOXELS_COUNT),
    geometry(NULL)//,
    //mpiCom(mpiCom_)
  { 
    //nLeaves=nVoxels=ROOT_VOXELS_COUNT;
    std::fill_n(nLeaves,MAX_TH,0);
    nUniqueVertices=0;
    nUniqueSegments=0;
    verticesAreAssigned=false;  
    segmentsAreAssigned=false;  
    initVoxels();
    initialize(x0,delta);
  }
  
  /** \brief default constructor.
   *  Bounding box needs to be set after construction by calling init(x0,delta)
   *  Default bounding box is the unit cube at origin
   */
  LocalAmrGridT(/*MpiCommunication *mpiCom_=glb::mpiComWorld*/):
    //rootVoxels(ROOT_VOXELS_COUNT),
    geometry(NULL)//,
    //mpiCom(mpiCom_)
  {  
    std::vector<double> x0(ND,0.0);
    std::vector<double> delta(ND,1.0);
    //nLeaves=nVoxels=ROOT_VOXELS_COUNT;
    std::fill_n(nLeaves,MAX_TH,0);
    nUniqueVertices=0;
    nUniqueSegments=0;
    verticesAreAssigned=false;
    segmentsAreAssigned=false;  
    initVoxels(); 
    initialize(&x0[0],&delta[0]);
  }

  ~LocalAmrGridT()
  {
    delete geometry;
  }
  
  /** \brief Reset the AMR grid to its initial state (i.e. a uniform grid at level 
   *  ROOT_LEVEL), deleting any refined voxel and returning their allocated memory 
   *  chunks to the memory pool.
   */
  void clear()
  {        
    //visitTree(internal::localAmrGridVisitor::ClearVoxelsT<MyType>(this),1);
    
    //#pragma omp parallel for num_threads(glb::num_omp_threads)
    for (long i=0;i<ROOT_VOXELS_COUNT;++i)
      rootVoxels[i].setEmpty();
    verticesAreAssigned=false;
    nUniqueVertices=0;
    nUniqueSegments=0;
    std::fill_n(nLeaves,MAX_TH,0);
    //nLeaves=nVoxels=ROOT_VOXELS_COUNT;
    for (int i=0;i<MAX_TH;++i)
      voxelGroupPool[i].freeChunks();
    //voxelGroupPool.freeChunks();
    segmentsAreAssigned=false;
  }

  /** \brief returns the number of unique vertices in the grid.
   *  \return the number of unique vertices if it is known, 0 otherwise.
   *  \note The number of unique vertices in the grid becomes known only after a call to 
   *  assignVerticesToLeaves() and it remains so until the grid is modified in any way
   *  (i.e. until voxels are recycled or created). 
   */ 
  unsigned long getUniqueVerticesCount()
  {
    if (verticesAreAssigned)
      return nUniqueVertices;
    else return 0;
  }

  /** \brief returns the number of unique segments in the grid.
   *  \return the number of unique segments if it is known, 0 otherwise.
   *  \note The number of unique segments in the grid becomes known only after a call to 
   *  assignVerticesToLeaves() and it remains so until the grid is modified in any way
   *  (i.e. until voxels are recycled or created).
   */ 
  unsigned long getUniqueSegmentsCount()
  {
    if (segmentsAreAssigned)
      return nUniqueSegments;
    else return 0;
  }

  /*
  void clear()
  {
  if (nVoxels!=ROOT_VOXELS_COUNT)
      {
	for (int i=0;i<ROOT_VOXELS_COUNT;++i)
	  {
	    clearRec(&rootVoxels[i]);    
	    rootVoxels[i].setEmpty();
	  }
	nLeaves=ROOT_VOXELS_COUNT;
	nVoxels=ROOT_VOXELS_COUNT;
      }
  }
  */  

  /** 
   *  \brief Assign each vertex (corner of voxels) and optionally segments to a single leaf
   *  \param assignSegments if true, also assign each segment to a single leaf
   *  \param nThreads number of threads to use 
   *  \param force If \a true, forces the function to reassign the vertices to voxels even
   *  in case the assignement has already been computed. By default (force=false), the 
   *  assignement is not recomputed if the structure has not been changed since last call 
   *  to assignVerticesToLeaves.
   *  \return the number of unique vertices in the grid
   *  \note The status of the vertices of a given voxel can then be queried
   *  using function Voxel::getVertexFlag() .
   */
  unsigned long assignVerticesToLeaves(bool assignSegments=false,
				       int nThreads=glb::num_omp_threads,
				       bool force=false)
  { 
    if ((segmentsAreAssigned||(!assignSegments))&&
	(verticesAreAssigned)&&
	(!force))
      return nUniqueVertices;
    
    if (nThreads<1) nThreads=glb::num_omp_threads;

    if (nThreads==1)
      {
	internal::localAmrGridVisitor::AssignVerticesT<MyType> 
	  visitor(this,assignSegments);

	visitTree(visitor,nThreads);
	nUniqueVertices = visitor.getNVertices();
	nUniqueSegments = visitor.getNSegments();
      }
    else
      {
	//internal::localAmrGridVisitor::AssignVerticesT<MyType> visitors[nThreads];
	std::vector< internal::localAmrGridVisitor::AssignVerticesT<MyType> >
	  visitors(nThreads);

	for (int i=0;i<nThreads;++i) visitors[i].init(this,assignSegments);
	visitTree(nThreads,&visitors[0]);
	nUniqueVertices=visitors[0].getNVertices();
	nUniqueSegments=visitors[0].getNSegments();
	for (int i=1;i<nThreads;++i)
	  {
	    nUniqueVertices+=visitors[i].getNVertices();
	    nUniqueSegments+=visitors[i].getNSegments();
	  }
      }
	
    glb::console->print<LOG_DEBUG>
      ("Assigned %ld unique vertices and %ld segments.\n",nUniqueVertices,nUniqueSegments);
	
    verticesAreAssigned = true;
    segmentsAreAssigned = assignSegments;


    return nUniqueVertices;
  }

  /** 
   *  \brief Returns a list of the segments a given voxel owns. Each segment may
   *  only be owned by a single voxel of the grid.
   *  \param v a pointer to the voxel to test
   *  \param out an outputIterator that can store std::pair<int,int> elements, where
   *   the first element is a vertex index and the second the index 'dim' of the dimension
   *   along which the segment is aligned. The vertices are always given such that 
   *   the segment goes in the positive direction along 'dim' 
   *  \return the number of owned segments
   *  \warning A call to assignVerticesToLeaves() with option assignSegments=true must be 
   *  performed before this function can return successfully
   */
  template <class OutputIterator>
  int getOwnedSegments(const Voxel *v, OutputIterator out) const
  {
    int n=0;
    unsigned int flags = v->getSegmentFlags();
  
    for (int i=0;i<Voxel::NVERT-1;++i)
      {	
	for (int j=0;j<NDIM;++j)
	  {
	    if (flags&1)
	      {
		(*out)=std::make_pair(i,j);
		++out;
		++n;
		//printf("Added one flag @(%d,%d) (->%d)\n",i,j,n);
	      }
	    flags>>=1;
	  }
      }
    
    return n;
  }

  /*
  template <class OutputIterator>
  int getOwnedVertices(Voxel *v, OutputIterator out) const
  {
    int n=0;    
    int flags=v->flags;
    for (int i=0;i<Voxel::NVERT;++i)
      {	
	if (flags&1)
	  {
	    (*out)=i;
	    ++out;
	    ++n;
	  }
	flags>>=1;
      }
    
    return n;
    }
  */ 

  /** \brief Builds the AMR grid from an unstructured mesh so that its local resolution is
   *   roughly the same as the mesh. The resolution is set such that any voxel intersecting 
   *   the bounding box of any given simplex in the mesh has a diagonal resolution smaller 
   *   than the smallest height of this simplex times the scale factor \p scaleFactor.
   *
   *  \tparam M  unstructured mesh type
   *  \param[in] mesh the unstructured mesh   
   *  \param     scaleFactor given a simplex, the size of the voxels it intersects is set 
   *             such that it is smaller than scaleFactor times smallest height of the 
   *             simplex.
   *  \param     maxLevel The maximum refinement level allowed
   *  \param     fromScratch rebuild the grid from scratch when true or try 
   *             to adapt the current grid when false. The second option will
   *             be faster if the current state was built from a similar mesh.
   * \param      nThreads Number of openMP threads to use
   * \warning when fromScratch=false, the amr grid will onyl be refined further where
   * needed, but never unrefeined.
   * \warning the simplices cache is erased
   */
  template <class M>
  void buildFromMesh(M *mesh, 
		     double scaleFactor=1.0,
		     char maxLevel=MAX_LEVEL_FROM_ROOT, 
		     bool fromScratch=true,
		     int nThreads=glb::num_omp_threads)
  {
    typedef typename M::vertexPtr_iterator vertexPtr_iterator;
    typedef typename M::ghostVertexPtr_iterator ghostVertexPtr_iterator;
    typedef typename M::simplexPtr_iterator simplexPtr_iterator;
    typedef typename M::ghostSimplexPtr_iterator ghostSimplexPtr_iterator;
    typedef typename M::simplexPtr_LG_iterator simplexPtr_LG_iterator;
    typedef typename M::Simplex Simplex;
    typedef typename M::GhostSimplex GhostSimplex;
    typedef typename M::Vertex Vertex;
    typedef typename M::SegmentHandle SegmentHandle;
    typedef typename M::Coord MCoord;
    typedef typename M::GeometricProperties MeshGeometricProperties;
    typedef typename Simplex::FacetHandle FacetHandle;

    //static const double NDIM_INV = 1.0/NDIM;
    //const double pixelDiagonalFactor = scaleFactor/sqrt(NDIM);
    const double pixelDiagonalFactor = scaleFactor;

    static bool refineRegionEdge=true; // must have LG iterators below if false

    typename TimerPool::Timer timer;
    timer.start();

    //char maxLevelTh[nThreads]={0};
    //char maxLevel=0;

    // FIXME: not implemented yet
    if (fromScratch) 
      clear();
    else
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("Option 'fromScratch=false' ignored (not implemented yet !)\n");
	clear();	
      }
        
    // Find the refinement level over each simplex, by computing the length
    // of it shortest height of each simplex
    /*
#pragma omp parallel for num_threads(nThreads)
    for (int k=0;k<nThreads;k++)
      {
	const simplexPtr_iterator its_end=mesh->simplexEnd();
	for (simplexPtr_iterator it=
	       mesh->simplexBegin(k,nThreads);it!=its_end;++it)
    */
    // This is better for load-balancing
    FOREACH_BATCH_SIMPLEX(mesh,nThreads,32,k,it)
      for (;it!=it_end;++it)
	{
	  Simplex *s=(*it);
	  // compute the simplex extent
	  double v=mesh->computeProjectedVolume(s);
	  double sMax=mesh->computeProjectedVolume(s->getFacetHandle(0).asPointer());
	  double sMin=sMax;

	  for (int i=1;i<Simplex::NFACET;++i)
	    {
	      double tmp=mesh->computeProjectedVolume(s->getFacetHandle(i).asPointer());
	      if (tmp>sMax) sMax=tmp;
	      if (tmp<sMin) sMin=tmp;
	    }
	    	    
	  double h=pixelDiagonalFactor * ((sMax==0)?0:(v*NDIM)/sMax);	 	  
	  char lvl=resolution2Level(h,maxLevel);
	 
	  // boundaries are refined at max level
	  if (M::BOUNDARY_TYPE != BoundaryType::PERIODIC)
	    {
	      for (int i=0;i<Simplex::NNEI;++i)
		{
		  if (s->getNeighbor(i)==NULL)
		    refineOverSimplex(mesh,s->getFacetHandle(i),maxLevel);
		}	      
	    }
	    
	  // MPI process boundaries are also refined at maxLevel
	  if (refineRegionEdge)
	    {
	      for (int i=0;i<Simplex::NNEI;++i)
		{
		  if ((s->getNeighbor(i)!=NULL)&&
		      (!s->getNeighbor(i)->isLocal()))
		    {
		      refineOverSimplex(mesh,s->getFacetHandle(i),maxLevel);
		    }
		}
	    }

	  refineOverSimplex(mesh,s,lvl);
	}
    FOREACH_BATCH_END
      //}

    // If we do not refine process boundaries at maxLevel, then we refine
    // the ghost simplices ...
    if (!refineRegionEdge)
      {
#pragma omp parallel for num_threads(nThreads)
	for (int k=0;k<nThreads;k++)
	  {
	    const ghostSimplexPtr_iterator itgs_end=mesh->ghostSimplexEnd();
	    for (ghostSimplexPtr_iterator it=
		   mesh->ghostSimplexBegin(k,glb::num_omp_threads);it!=itgs_end;++it)
	      {
		GhostSimplex *s=*it;	  
		
		double v=mesh->computeProjectedVolume(s);
		double sMax=mesh->computeProjectedVolume(s->getFacetHandle(0).asPointer());
		double sMin=sMax;

		for (int i=1;i<Simplex::NFACET;++i)
		  {
		    double tmp=
		      mesh->computeProjectedVolume(s->getFacetHandle(i).asPointer());

		    if (tmp>sMax) sMax=tmp;
		    if (tmp<sMin) sMin=tmp;
		  }
	    	    
		double h=pixelDiagonalFactor*((sMax==0)?0:(v*NDIM)/sMax);	    	    
		char lvl=resolution2Level(h,maxLevel);

		if (M::BOUNDARY_TYPE != BoundaryType::PERIODIC)
		  {
		    for (int i=0;i<Simplex::NNEI;++i)
		      {
			if (s->getNeighbor(i)==NULL)
			  refineOverSimplex(mesh,s->getFacetHandle(i),maxLevel);
		      }	      
		  }
		refineOverSimplex(mesh,s,lvl);
	      }
	  }
      }
  }  

private:
  class DummySimplexTagger{
      template <class AMR, class Simplex>
      void tagSimplex(AMR *amr,Simplex *s, 
		      double bBox[2][Simplex::NDIM],double maxLevelVoxelLengthInv[Simplex::NDIM],
		      double v, double sMin, double sMax) const
      {}
  };

public:
  /** \brief Builds the AMR grid from an unstructured mesh so that its local resolution is
   *   roughly the same as the mesh. This version is faster but its root voxel can only
   *   be refine at a constant level (level may differ between roots, but not within 
   *   each root)
   *
   * \tparam M  unstructured mesh type
   * \param[in] mesh the unstructured mesh    
   * \param     scaleFactor given a simplex, the size of the voxels it intersects is set 
   *             such that it is smaller than scaleFactor times smallest height of the 
   *             simplex
   * \param     maxLevel The maximum refinement level allowed
   * \param     fromScratch rebuild the grid from scratch when true or try 
   *             to adapt the current grid when false. The second option will
   *             be faster if the current state was built from a similar mesh.
   * \param      nThreads Number of openMP threads to use
   * \param      assignVertices If true assign vertices uniquely to voxels (faster 
   *             than calling assignVerticesToLeaves later)
   * \param      assignSegments If true assign segments uniquely to voxels (faster 
   *             than calling assignVerticesToLeaves later)  
   *             
   * \warning when fromScratch=false, the amr grid will only be refined further where
   * needed, but never unrefeined.
   * \warning the simplices cache is erased
   */
  template <class M>
  void buildFromMesh_Fast(M *mesh, 
			  double scaleFactor=1.0,
			  char maxLevel=MAX_LEVEL_FROM_ROOT, 			  
			  bool fromScratch=true,
			  int nThreads=glb::num_omp_threads,
			  bool assignVertices=false,
			  bool assignSegments=false)
  {  
    DummySimplexTagger dummy;
    buildFromMesh_Fast(mesh,dummy,scaleFactor,
		       maxLevel,fromScratch,nThreads,
		       assignVertices,assignSegments);
  }


  /** \brief Builds the AMR grid from an unstructured mesh so that its local resolution is
   *   roughly the same as the mesh. This version is faster but its root voxel can only
   *   be refine at a constant level (level may differ between roots, but not within 
   *   each root)
   *
   * \tparam M  unstructured mesh type
   * \param[in] mesh the unstructured mesh  
   * \param[in] tagger a functor that can tag mesh simplices as they are scanned according
   * to their configuration (used to set special cases in LocalAmrGridProjector)
   * \param     scaleFactor given a simplex, the size of the voxels it intersects is set 
   *             such that it is smaller than scaleFactor times smallest height of the 
   *             simplex
   * \param     maxLevel The maximum refinement level allowed
   * \param     fromScratch rebuild the grid from scratch when true or try 
   *             to adapt the current grid when false. The second option will
   *             be faster if the current state was built from a similar mesh.
   * \param      nThreads Number of openMP threads to use
   * \param      assignVertices If true assign vertices uniquely to voxels (faster 
   *             than calling assignVerticesToLeaves later)
   * \param      assignSegments If true assign segments uniquely to voxels (faster 
   *             than calling assignVerticesToLeaves later)
   *             
   * \warning when fromScratch=false, the amr grid will only be refined further where
   * needed, but never unrefeined.
   * \warning the simplices cache is erased
   */
  template <class M, class SimplexTagger>
  void buildFromMesh_Fast(M *mesh, const SimplexTagger &tagger,
			  double scaleFactor=1.0,
			  char maxLevel=MAX_LEVEL_FROM_ROOT, 			  
			  bool fromScratch=true,
			  int nThreads=glb::num_omp_threads,
			  bool assignVertices=false,
			  bool assignSegments=false)
  //bool tagSimplicesForProjectionsetSimplexCache=true)
  {
    typedef typename M::vertexPtr_iterator vertexPtr_iterator;
    typedef typename M::ghostVertexPtr_iterator ghostVertexPtr_iterator;
    typedef typename M::simplexPtr_iterator simplexPtr_iterator;
    typedef typename M::ghostSimplexPtr_iterator ghostSimplexPtr_iterator;
    typedef typename M::simplexPtr_LG_iterator simplexPtr_LG_iterator;
    typedef typename M::Simplex Simplex;
    typedef typename M::GhostSimplex GhostSimplex;
    typedef typename M::Vertex Vertex;
    typedef typename M::SegmentHandle SegmentHandle;
    typedef typename M::Coord MCoord;
    typedef typename M::GeometricProperties MeshGeometricProperties;
    typedef typename Simplex::FacetHandle FacetHandle;
    
    //static const double NDIM_INV = 1.0/NDIM;
    //const double pixelDiagonalFactor = scaleFactor/sqrt(NDIM);
    const double pixelDiagonalFactor = scaleFactor;
    static bool refineRegionEdge=true; // must have LG iterators below if false

    typename TimerPool::Timer timer;
    timer.start();
    
    // FIXME: not implemented yet
    if (fromScratch) 
      clear();   
    else
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>
	  ("Option 'fromScratch=false' ignored (not implemented yet !)\n");
	clear();
	//exit(-1);
      }
    
    double minVoxelVolume=1.0;
    double maxLevelVoxelLengthInv[NDIM];    
    for (int i=0;i<NDIM;++i)
      {
	maxLevelVoxelLengthInv[i]=(1.0/getVoxelLength(maxLevel,i));
	minVoxelVolume*=getVoxelLength(maxLevel,i);
      }
    
    double maxLevelVoxelLength_max=getVoxelLength(maxLevel,0);
    for (int i=1;i<NDIM;++i)
      if (maxLevelVoxelLength_max<getVoxelLength(maxLevel,i))
	maxLevelVoxelLength_max=getVoxelLength(maxLevel,i);

    //double simplexVolumeThreshold=minVoxelVolume*pow(1.e-3,NDIM);
    //double simplexAnisotropyThreshold=1<<(ROOT_LEVEL+maxLevel+3);   

    /* 
    double volumeRatio = pow(100,-NDIM);
    double volumeThreshold=getVoxelVolume(maxLevel)*volumeRatio;
    double surfaceThreshold=pow(getVoxelVolume(maxLevel)*volumeRatio,2.0/NDIM);
    double lengthThreshold=pow(getVoxelVolume(maxLevel)*volumeRatio,1.0/NDIM);
    */
    /*
    double refLevel = std::max((int)maxLevel,(int)(10-ROOT_LEVEL));
    double volumeRatio = pow(1000,-NDIM);
    double volumeThreshold=getVoxelVolume(refLevel)*volumeRatio;
    //double surfaceThreshold=pow(getVoxelVolume(maxLevel)*volumeRatio,2.0/NDIM);
    double lengthThreshold=pow(getVoxelVolume(refLevel)*volumeRatio,1.0/NDIM);
    double anisotropyThreshold=100;
    */
    std::vector< std::vector<unsigned char> > rootLevel(nThreads);
    for (int i=0;i<nThreads;++i)
      rootLevel[i].assign(ROOT_VOXELS_COUNT,0);

    bool setSimplexCache = !hlp::SameType<SimplexTagger,DummySimplexTagger>::value;

    if ((maxLevel>0)||(setSimplexCache))
      {
	// This is better for load-balancing
	FOREACH_BATCH_SIMPLEX(mesh,nThreads,32,k,it)
	  for (;it!=it_end;++it)
	    {
	      int th=omp_get_thread_num();
	      Simplex *s=(*it);
	      // compute the simplex extent
	      double v=mesh->computeProjectedVolume(s);
	      double sMax=mesh->computeProjectedVolume(s->getFacetHandle(0).asPointer());
	      double sMin=sMax;
	      //double anisotropy=sMax/sMin;
	      //double hMin=v/sMax;

	      for (int i=1;i<Simplex::NFACET;++i)
		{
		  double tmp=mesh->computeProjectedVolume(s->getFacetHandle(i).asPointer());
		  if (tmp>sMax) sMax=tmp;
		  if (tmp<sMin) sMin=tmp;
		}
	    	    
	      double h=pixelDiagonalFactor*((sMax==0)?0:(v*NDIM)/sMax);
	      char reqLvl=resolution2Level(h,maxLevel+1);
	      char lvl = std::min(reqLvl,maxLevel);

	      double bBox[2][NDIM];
	      mesh->template computeBoundingBox<Simplex*,NDIM>(s,bBox);	          
	     	      
	      tagger.tagSimplex(this,s,bBox,maxLevelVoxelLengthInv,v,sMin,sMax);
	      /*
	      if (setSimplexCache)
		{	
		  s->cache.c[0]=ProjectionTag::regular;
		  // If this is true, simplex is smaller than half the voxels size, so
		  // check if it falls entirely within a single voxel. If So, tag it by
		  // setting cache.c[] to 1;
		  if (reqLvl>=maxLevel)
		    {
		      bool inside=true;//ProjectionTag::sampleNoOverlap;//(1<<0);
		      for (int j=0;j<NDIM;++j)
			{		  
			  int a=(int)((bBox[0][j]-getBBoxMin(j)-getEpsilon(j))*
				      maxLevelVoxelLengthInv[j]);
			  int b=(int)((bBox[1][j]-getBBoxMin(j)+getEpsilon(j))*
				      maxLevelVoxelLengthInv[j]);
			  inside &= (a==b);
			}
		      if (inside) s->cache.c[0]=ProjectionTag::sampleNoOverlap;//inside;
		      //printf("Tagged !\n");
		    }
		  
		  if ((sMax > sMin*anisotropyThreshold)|| //anisotropy
		      (v<volumeThreshold)|| //volume
		      (v<sMax*lengthThreshold)) // minimum width
		    {
		      s->cache.c[0]|=ProjectionTag::sampleWithOverlap; // (1<<1)
		    }
	
		}
	      */
	      if (maxLevel==0) continue;
		
	      internal::localAmrGridOverlapVisitor::SetRootLevelT<MyType> 
		visitor(this,lvl,&rootVoxels[0],rootLevel[th]);		  
	      visitBBoxOverlap(bBox,visitor);
	
	      visitor.setLevel(maxLevel);
	      // boundaries are refined at max level
	      if (M::BOUNDARY_TYPE != BoundaryType::PERIODIC)
		{
		  for (int i=0;i<Simplex::NNEI;++i)
		    {
		      if (s->getNeighbor(i)==NULL)
			{
			  mesh->template computeBoundingBox<FacetHandle,NDIM>
			    (s->getFacetHandle(i),bBox);
			  visitBBoxOverlap(bBox,visitor);
			  //refineOverSimplex(mesh,s->getFacetHandle(i),maxLevel);
			}
		    }	      
		}
	    
	      // MPI process boundaries are also refined at maxLevel
	      if (refineRegionEdge)
		{
		  for (int i=0;i<Simplex::NNEI;++i)
		    {
		      if ((s->getNeighbor(i)!=NULL)&&
			  (!s->getNeighbor(i)->isLocal()))
			{		     
			  mesh->template computeBoundingBox<FacetHandle,NDIM>
			    (s->getFacetHandle(i),bBox);
			  visitBBoxOverlap(bBox,visitor);
			  //refineOverSimplex(mesh,s->getFacetHandle(i),maxLevel);
			}
		    }
		}
	    }
	FOREACH_BATCH_END;
     

#pragma omp parallel for num_threads(nThreads) schedule(dynamic,1)
	for (long i=0;i<ROOT_VOXELS_COUNT;++i)
	  {
	    unsigned char lvl=rootLevel[0][i];
	    for (int th=1;th<nThreads;++th)
	      if (lvl < rootLevel[th][i]) lvl=rootLevel[th][i];
	
	    rootVoxels[i].refine(this,lvl);
	    rootLevel[0][i]=lvl;
	  }      
      }

    if (assignVertices||assignSegments)
      {
	if (nThreads==1)
	  {
	    internal::localAmrGridVisitor::AssignVertices_FastT<MyType> 
	      visitor(this,assignSegments,rootVoxels,&rootLevel[0][0]);

	    visitTree(visitor,nThreads);
	    nUniqueVertices = visitor.getNVertices();
	    nUniqueSegments = visitor.getNSegments();
	  }
	else
	  {
	    //internal::localAmrGridVisitor::AssignVerticesT<MyType> visitors[nThreads];
	    std::vector< internal::localAmrGridVisitor::AssignVertices_FastT<MyType> >
	      visitors(nThreads);

	    for (int i=0;i<nThreads;++i) 
	      visitors[i].init(this,assignSegments,rootVoxels,&rootLevel[0][0]);

	    visitTree(nThreads,&visitors[0]);
	    nUniqueVertices=visitors[0].getNVertices();
	    nUniqueSegments=visitors[0].getNSegments();
	    for (int i=1;i<nThreads;++i)
	      {
		nUniqueVertices+=visitors[i].getNVertices();
		nUniqueSegments+=visitors[i].getNSegments();
	      }
	  }

	glb::console->print<LOG_DEBUG>
	  ("Assigned %ld unique vertices and %ld segments.\n",nUniqueVertices,nUniqueSegments);
	verticesAreAssigned = true;
	segmentsAreAssigned = assignSegments;
      }
    /*
     FOREACH_BATCH_SIMPLEX(mesh,nThreads,32,k,it)
      for (;it!=it_end;++it)
	{
	  Simplex *s=(*it);
	  if (s->cache.c[0])
	    {
	      Data *data=&(getVoxelAt(s->getVertex(0)->getCoordsPtr())->data);
	      Data val=s->mass.getValue();
#pragma omp atomic
	      (*data)+=val;
	    }
	}
    FOREACH_BATCH_END;
    */
  }  


  /** \brief Converts a resolution to corresponding refinement level in the AMR grid. 
   *
   *  The refinement level is computed so that a voxel at this level has 
   *  size S such that resolution/2 < S < resolution if possible.
   *  \param resolution the required resolution 
   *  \param maxLevel the maximum refinement level allowed
   *  \return the level corresponding to resolution or maxLevel if maxLevel is smaller 
   */
  int resolution2Level(double resolution, int maxLevel) const
  {
    ICoord ratio = deltaX_max/resolution;        
    int level;

    // Compute the level of refinement corresponding to 'resolution'    
    if (ratio >= (1L<<ROOT_LEVEL<<maxLevel)) 
      {
	level=maxLevel;
      }
    else
      {
	ratio >>=ROOT_LEVEL;
	for (level=0;(ratio>>level);level++){}
      }   

    return level;
  }

  /** \brief Converts a resolution to corresponding refinement level in the AMR grid. 
   *  The refinement level is computed so that a voxel at this level has 
   *  size S such that resolution/2 < S < resolution if possible.
   *  \param resolution the required resolution 
   *  \return the level corresponding to resolution 
   */
  int resolution2Level(double resolution) const
  {
    ICoord ratio = deltaX_max/resolution;        
    int level;

    ratio >>=ROOT_LEVEL;
    for (level=0;(ratio>>level);level++){}
      
    return level;   
  }

  /** \brief Refine the AMR grid if needed to ensure that, if possible, the size along
   *   each axis of the leaf containing the point at coordinates 'coords' is at 
   *   least S < resolution.
   *  \param coords coordinates of the point
   *  \param level The refinement level needed at point \a coord
   *  \param maxLevel the maximum refinement level allowed
   */
  template <class InputIterator>
  Voxel *setMinimumLevelAt(InputIterator coords, int level, 
			   int maxLevel=MAX_LEVEL_FROM_ROOT)
  {           
    ICoord index = coords2Index(coords);
    Voxel *voxel = getVoxel(index,level);      
    int nRef = level-voxel->getLevel();

    for (int i=0;i<nRef;++i)
      {
	refine(voxel);
	voxel=voxel->getChild(voxel->getQuadrant(index));
      }
   
    return voxel;
  }
  
  /** \brief Projects the weight function defined over the unstructured mesh \a mesh 
   *  onto the AMR grid. The Function is defined at order 0 (i.e. constant over simplices)
   * \param[in] mesh the unstructured mesh
   * \param[in] wf a functor that returns the weight of a SIMPLEX given as 
   * argument (see WF)  
   * \param[in] checkTags if true, read the simplices cache.c[0] value and interpret it as 
   * follows: 0: nothing special, 1: simplex is fully contained inside a voxel
   * 2 or 3: simplex is very close to degenerate
   * \param[in] nThreads the number of threads to use, -1 to use as many as possible   
   * \param[in] verbose if true, progression and timings will be printed to the console
   * \return the number of simplices that had to be reprojected
   * \tparam M  unstructured mesh type
   * \tparam WF the class of the weight functor that implements a 
   * WF::operator()(MESH::Simplex *s) method returning the weight of a simplex \a s.
   * \tparam checkAccuracy Enable/Disable accuracy checking
   * \tparam IF The floating point type to use for intermediate computation. Defaults to
   * long double, see LocalAmrGridProjectorT.
   * \tparam HF   A high precision floating point type to use for the computation of 
   * the individual contributions to each voxel when IF is not precise enough. Defaults 
   * to 'long double'.
   * \tparam SF   A floating point type to use for the summation of all the individual 
   * contributions to each voxel. If this is different from Amr::Data, then an additional 
   * temporary array of size sizeof(SF)*number_of_voxels will have to be allocated. 
   * Defaults to 'long double'.       
   */
  template <class M,class WF, bool checkAccuracy,
	    class IF=long double,
	    class HF=long double,
	    class SF=long double>
  long projectMesh0(M *mesh, const WF &wf,
		    double accLevel=checkAccuracy?0.1:0,
		    bool checkTags=true,
		    int nThreads=glb::num_omp_threads,
		    bool verbose=false)
  {
    typedef LocalAmrGridProjectorT<MyType,M,checkAccuracy,IF,HF,SF> Projector;
    Projector projector(this,mesh,accLevel,nThreads,verbose);
    return projector.template project<WF>(wf,checkTags);   
  }

  /** \brief Projects the weight function defined over the unstructured mesh \a mesh 
   *  onto the AMR grid. The Function is defined at order 1 (i.e. linearly interpolated 
   *  over the simplices)   
   * \param[in] mesh the unstructured mesh
   * \param[in] wf a functor that returns the value of the weight at VERTICES locations 
   *  (vertex is given as argument to the functor, which returns its weight, see WF)
   * \param[in] wdf a functor that sets the derivative of the weight function over
   * each SIMPLEX given as argument (see WDF).
   * \tparam checkAccuracy Enable/Disable accuracy checking
   * \param[in] checkTags if true, read the simplices cache.c[0] value and interpret it as 
   * follows: 0: nothing special, 1: simplex is fully contained inside a voxel
   * 2 or 3: simplex is very close to degenerate
   * \param[in] nThreads the number of threads to use, -1 to use as many as possible   
   * \param[in] verbose if true, progression and timings will be printed to the console
   * \return the number of simplices that had to be reprojected
   * \tparam M  unstructured mesh type
   * \tparam WF The class of the weight functor that implements a 
   * WF::operator()(MESH::Vertex *v) method returning the weight evaluated at the location
   * of vertex \a v.   
   * \tparam WDF the class of the weight derivative functor that implements a 
   * WDF::operator()(MESH::Simplex *s, double result[NDIM]) method that sets the NDIM values
   * of result to the gradient of the weight within the simplex \a s. The gradient
   * should be computed such that WF(V_i) = WF(V_0) + grad(W) * (V_i-V_0), whith V_i the
   * ith vertex of the current simplex.
   * \tparam IF The floating point type to use for intermediate computation. Defaults to
   * long double, see LocalAmrGridProjectorT
   * \tparam HF   A high precision floating point type to use for the computation of 
   * the individual contributions to each voxel when IF is not precise enough. Defaults 
   * to 'long double'.
   * \tparam SF   A floating point type to use for the summation of all the individual 
   * contributions to each voxel. If this is different from Amr::Data, then an additional 
   * temporary array of size sizeof(SF)*number_of_voxels will have to be allocated. 
   * Defaults to 'long double'.    
   */
  template <class M, class WF, class WDF, bool checkAccuracy,
	    class IF=long double,
	    class HF=long double,
	    class SF=long double>
  long projectMesh1(M *mesh,const WF &wf,const WDF &wdf,
		    double accLevel=checkAccuracy?0.1:0,
		    bool checkTags=true,
		    int nThreads=glb::num_omp_threads, 
		    bool verbose=false)
  {
    typedef LocalAmrGridProjectorT<MyType,M,checkAccuracy,IF,HF,SF> Projector;
    Projector projector(this,mesh,accLevel,nThreads,verbose);
    return projector.template project<WF,WDF>(wf,wdf,checkTags);   
  }
  
  /** \brief Retrieve a voxel with index 'index' if its level is less or equal to maxLevel
   *  or its parent at level maxLevel. 
   *  \param index the index of the voxel
   *  \param maxLevel The maximum level the retrieved voxel can have
   *  \return the voxel with index 'index' or it parent at level maxLevel
   */
  Voxel *getVoxel(ICoord index, int maxLevel = MAX_LEVEL_FROM_ROOT)
  {
    ICoord rootIndex=getICoordFromIndex(index,0,0);
    int dec=ROOT_LEVEL;
    for (int i=1;i<NDIM;++i)
      {
	rootIndex += getICoordFromIndex(index,i,0)<<dec;
	dec+=ROOT_LEVEL;
      }
   
    Voxel *voxel = &rootVoxels[rootIndex];
    
    while ((!voxel->isLeaf())&&
	   (voxel->getLevel()<maxLevel)&&
	   (index!=voxel->getIndex()))
      {
	/*
	if (getQuadrantAtLevel(index, voxel->getLevel()) != voxel->getQuadrant(index))
	  printf("Quadrant differ : %d / %d \n",
		 getQuadrantAtLevel(index, voxel->getLevel()),
		 voxel->getQuadrant(index));
	*/
	voxel = voxel->getChild(voxel->getQuadrant(index));	
      }
          
    return voxel;
  }
  /*
  template <class TT>
  bool CheckInBound(const TT * x)
  {
    if (geometry->onBoundary())
      {
	
      }
  }
  */
  
  /** \brief Retrieve the most refined voxel in the grid with level at most maxLevel
   *  that contains the points with coordinates 'coords'.
   *  \param coords the coordinates of the point
   *  \param includeRightBoundary if true, then return NULL also if the point is exactly 
   *  on the right boundary (i.e. according to simulation of simplicity convention, the
   *  rightmost boundary should NOT be within the box for non periodic boundary conditions)
   *  \param maxLevel the maximum level the retieved voxel can reach
   *  \return the leaf voxel containing point coords or its parent at level maxLevel. If 
   *  coords fall exactly on the boundary of several voxels, the returned voxel is the
   *  one for which the coordinates of its center are each maximal. If coords are out of
   *  bound and boudnary conditions are non periodic, NULL is returned.
   *  \see getVoxelAt_noCheck
  */
  template <class InputIterator>
  Voxel *getVoxelAt(InputIterator coords, bool includeRightBoundary=true, 
		    int maxLevel = MAX_LEVEL_FROM_ROOT)
  {
    int warning;
    ICoord index = coords2Index(coords,warning,includeRightBoundary);
    if (warning) return NULL;
    return getVoxel(index,maxLevel);    
  }

  /** \brief Retrieve the most refined voxel in the grid with level at most maxLevel
   *  that contains the points with coordinates 'coords'. No test on the validity of coords
   *  is performed.
   *  \param coords the coordinates of the point
   *  \param maxLevel the maximum level the retieved voxel can reach
   *  \return the leaf voxel containing point coords or its parent at level maxLevel. If 
   *  coords fall exactly on the boundary of several voxels, the returned voxel is the
   *  one for which the coordinates of its center are each maximal.
   *  \see getVoxelAt
  */
  template <class InputIterator>
  Voxel *getVoxelAt_noCheck(InputIterator coords, int maxLevel = MAX_LEVEL_FROM_ROOT)
  {
    ICoord index = coords2Index_noCheck(coords);
    return getVoxel(index,maxLevel);    
  }
  
  /**
   * \brief returns an integer corresponding to the quadrant (octant in 3D) a voxel with 
   * index 'voxelIndex' falls compared to the center of a refernce voxel with index 'ref'.
   * \param ref  the index of the reference voxel
   * \param voxelIndex the index of the voxel we are trying to locate
   * \return an integer Q ranging from 0 to (1<<NDIM). If voxel with index 'voxelIndex' has 
   * a higher coordinate along dimension D then the Dth bit of Q is 1, and 0 otherwise.  
   */  
  static int getQuadrantFromRef(ICoord ref, ICoord voxelIndex)
  {    
    int result=0;    
    if ((ref&INDEX_MASK)<(voxelIndex&INDEX_MASK)) result ++;
   
    //if (ref!=voxelIndex)
      {
	for (int i=1;i<NDIM;i++)
	  {
	    ref   >>= INDEX_DEC;
	    voxelIndex >>= INDEX_DEC;
	
	    if ((ref&INDEX_MASK)<(voxelIndex&INDEX_MASK)) result += (1L<<i);	    
	  }
      }
      //else result=-1;

    return result;
  }
  

  /**
   * \brief returns an integer corresponding to the quadrant (octant in 3D) within its 
   * parent cell a voxel with index 'voxelIndex' and level 'voxelIndexLevel' would fall in.
   * This is faster than getQuadrantFromRef.
   * \param voxelIndex the index of the reference voxel
   * \param voxelIndexLevel the level of the voxel with index 'voxelIndex'
   * \return an integer Q ranging from 0 to (1<<NDIM) representing in which quadrant of its
   * parent cell the voxel with index 'voxelIndex' falls. The kth bit of Q is set to 0/1
   * if it lies on the left/right of its parent cell along dimension k.
   */
  static int getQuadrantAtLevel(ICoord voxelIndex, int voxelIndexLevel)
  {        
    ICoord dec =(MAX_LEVEL_FROM_ROOT - voxelIndexLevel);
    const ICoord testBit = 1 << dec;
    int result= ((voxelIndex&INDEX_MASK)&testBit)>>dec;
    
    for (int i=1;i<NDIM;++i)
      {
	dec--;
	voxelIndex >>= INDEX_DEC;
	result |= ((voxelIndex&INDEX_MASK)&testBit)>>dec;	    
      }
  
    return result;
  }

  /** \brief Refine a voxel, creating (1<<NDIM) child.
   *  \param voxel a pointer to the voxel to refine
   *  \return A pointer to the voxel that was refined (='voxel')
   */
  Voxel *refine(Voxel *voxel)
  {
    return voxel->refine(this);
  }

  /** \brief Sets the grid's bounding box. This needs to be called if the class was 
   *  initialized with the default constructor and a specific bounding box other
   *  than the unit cube with corner at origin needs to be used
   *  \param x0 an iterator to the lower left coordinates of the bounding box
   *  \param delta an iterator to the size of the bounding box along each dimension
   *  \warning Simulation of simplicity will only be guaranteed to work if the box origin and size can
   *  be expressed exactly, so it is usualy a good idea to use integral values for x0 and delta or numbers
   *  that can be expressed as sums of inverse powers of 2 (e.g 0.5 is fine but 0.2 is risky ...)
   */
  template <class InputIterator>
  void initialize(InputIterator x0, InputIterator delta)//, MpiCommunication *mpiCom_=NULL)
  {  
    //if (mpiCom_ != NULL) mpiCom = mpiCom_;
    xMin.assign(x0,x0+ND);
    deltaX.assign(delta,delta+ND);
    xMax.resize(ND);
    deltaX_inv.resize(ND);
    ICoordUnitLen.resize(ND);
    epsilon.resize(ND);
    halfEpsilon.resize(ND);
    halfEpsilonILen.resize(ND);
    deltaX_max=deltaX_min=deltaX[0];
    rootIStride.assign(ND+1,1);
    for (int i=0;i<ND;++i)
      {
	xMax[i] = xMin[i]+deltaX[i];
	deltaX_inv[i]=1./deltaX[i];
	if (deltaX[i]>deltaX_max) deltaX_max=deltaX[i];
	if (deltaX[i]<deltaX_min) deltaX_min=deltaX[i];
	ICoordUnitLen[i]=deltaX[i]/BBOX_ILEN;
	//printf("ICoordUnitLen[%d]=%g\n",i,ICoordUnitLen[i]);
	epsilon[i]=std::max(fabs(xMin[i]),fabs(xMax[i]))*
	  (std::numeric_limits<double>::epsilon()*20.0);

	//epsilon[i]=deltaX[i]*(std::numeric_limits<double>::epsilon()*100.0);
	halfEpsilon[i]=epsilon[i]/2;
	halfEpsilonILen[i]=halfEpsilon[i]*BBOX_ILEN;
	//ICoordUnitLen[i]/4;
	rootIStride[i+1]=N_ROOT_PER_DIM*rootIStride[i];
      }   

    for (int i=0;i<=MAX_LEVEL_FROM_ROOT;++i)
      {
	voxelVolume[i]=getVoxelLength(i,0);
	for (int j=1;j<NDIM;++j)
	  voxelVolume[i]*=getVoxelLength(i,j);

	voxelInverseVolume[i]=1.0/voxelVolume[i];
      }

#ifdef HAVE_BOOST
#ifdef HAVE_GMP
    // multiprecision variable initialization
    mp_xMin.assign(x0,x0+ND);
    mp_deltaX.assign(delta,delta+ND);
    mp_xMax.resize(ND);
    mp_deltaX_inv.resize(ND);
    mp_BBOX_ILEN = +BBOX_ILEN;
    for (int i=0;i<ND;++i)
      {
	mp_xMax[i]=mp_xMin[i]+mp_deltaX[i];
	mp_deltaX_inv[i]=1.0;
	mp_deltaX_inv[i] /= mp_deltaX[i];
      }
#endif
#endif
    
    clear();
    delete geometry;
    geometry = new GeometricProperties(x0,delta);
  }
  
  /** \brief return the maximum size_max of the bounding box along any direction.
   *  \return value of size_max
   */
  double getMaxLength() const
  {
    return deltaX_max;
  }

  /** \brief return a value epsilon such that epsilon is smaller that half the size of any
   *   voxel along dimension dim.
   *  \param dim the dimension along which to compute epsilon
   *  \return the value of epsilon
   */
  double getEpsilon(int dim) const
  {
    return epsilon[dim];
  }

  /** \brief return the volume of a voxel at level \a level
   *  \param level the level of the voxel
   *  \return the volume of the voxel
   */
  double getVoxelVolume(int level) const
  {
    return voxelVolume[level];
  }

  /** \brief return the inverse of the volume of a voxel at level \a level
   *  \param level the level of the voxel
   *  \return the inevrse of the volume of the voxel
   */
  double getVoxelInverseVolume(int level) const
  {
    return voxelInverseVolume[level];
  }
 
  /** \brief get the number of leaves in the grid.
   *  \return the number of leaves. This number is always greater than or equal to
   *  ROOT_VOXELS_COUNT (i.e. the number of root leaves that cannot be unrefined)
   * \warning This takes MAX_TH additions, so cache it to use it in a for loop ;)
   */
  unsigned long getNLeaves() const
  {
    long result=ROOT_VOXELS_COUNT;
    for (int i=0;i<MAX_TH;++i)
      result+=nLeaves[i];
    return result;
  }
  

  /** \brief Retrieve the number of voxels at a given \a level.
   *  \param level the level of the voxels or a negative number
   *  \return the number of voxels allocated at level \a level or the total number
   *  of voxels in the grid if \a level<0. 
   * \warning NOT IMPLEMENTED !
   */
  unsigned long getNVoxels(int level) const;
  /*
  unsigned long getNVoxels(int level=-1) const
  {
    // NOT IMPLEMENTED
    
    unsigned long result=ROOT_VOXELS_COUNT;
    if (level<0) result=nVoxels;
    else 
      {
	result=voxelGroupPool[level].getUsedCount()<<NDIM;
	//result=voxelGroupPool.getUsedCount()<<NDIM;
      }
    
    return result;
    
  }
*/
  /** \brief Retrieve the number of voxels
   *  \return the total number of voxels allocated (including the ROOT_VOXELS_COUNT 
   *  root voxels).
   * \warning This takes MAX_TH additions, so cache it to use it in a for loop ;)
   */
  unsigned long getNVoxels() const
  {
    unsigned long result=ROOT_VOXELS_COUNT;
    for (int i=0;i<MAX_TH;++i)
      result=voxelGroupPool[i].getUsedCount()<<NDIM;
   
    return result;
  }

  /** \brief Save the AMR grid to an ASCII file
   *  \param fname The name of the file, or NULL for default name
   */
  void toAscii(const char *fname=NULL)
  {
    static int count =0;
    char name[255];

    if (fname==NULL)
      sprintf(name,"amr_%6.6d.dat",count++);
    else 
      sprintf(name,"%s.amr.ascii",fname);

    FILE *f=fopen(name,"w");
    
    visitTree(internal::localAmrGridVisitor::ToAsciiT<MyType>(this,f),1);
    /*
    for (int i=0;i<ROOT_VOXELS_COUNT;++i)
      toAsciiRec(&rootVoxels[i],f);
    */
    fclose(f);
  }  

  /** \brief Save the AMR grid to a VTK file
   *  \param fname The name of the file, or NULL for default name
   */
  void toVtk(const char *fname=NULL)
  {
    IO::VtkAmrWriterT<MyType,float> vtkWriter(this,fname);
    vtkWriter.write();
  }

  /** \brief returns the geometry corresponding to the bounding box
   *  \return a pointer to the geometry class of the AMR grid
   */
  GeometricProperties *getGeometry()
  {
    return geometry;
  }

  /** \brief Get the orientation of a given direction along dimension 'dim'. For instance,
   *  getOrientation(Left,0)=-1, getOrientation(Left,1)=0, getOrientation(Bottom,2)=-1 and 
   *  getOrientation(Front|Top|Left,2)=1.
   *  \param dir the direction
   *  \param dim which dimension to check
   *  \return {-1,0,+1} if the direction is negative / undetermined / positive 
   *  along dimension 'dim'   
   */
  static int getOrientation(Direction dir, int dim)
  {
    return ((dir&(1L<<(16+dim)))>>(15+dim))-1;
  }

  /** \brief Get the direct neighbor of voxel 'v' at the same level or its highest 
   *  level parent if it does not exist.
   *  \param v A pointer to the voxel
   *  \param dir The direction toward the neighbor
   *  \return the neighbor voxel (or one of its parent); or 'v' itself if the neighbor is
   *  not within the non periodic bounding box.
   *  \warning untested
   */  
  Voxel *getNeighbor(Voxel *v, Direction dir) const
  {
    ICoord index = getNeighborIndex(v->getIndex(),v->getLevel(),dir);
    return getVoxel(index,v->getLevel());
  }

  /** \brief Get the direct neighbor of voxel 'v' at the same level or its highest 
   *  level parent if it does not exist. The boundary conditions are considered
   *  to be periodic whatever the configuration of the AMR grid.
   *  \param v A pointer to the voxel
   *  \param dir The direction toward the neighbor
   *  \return the neighbor voxel (or one of its parent)
   *  \warning untested
   */  
  Voxel *getNeighbor_periodic(Voxel *v, Direction dir) const
  {
    ICoord index = getNeighborIndex_periodic(v->getIndex(),v->getLevel(),dir);
    return getVoxel(index,v->getLevel());
  }

  /** \brief Get the index of the same level direct neighbor of a voxel with index 'index'
   *  and level 'level'. Note that these voxel do not have to exist.
   *  \param index the index whose neighbor we look for 
   *  \param level the level of the current index
   *  \param dir the direction toward the neighbor
   *  \return the index of the neighbor the parameter 'index' if the neighbor is not 
   *  within the bounding box and boundary conditions are not periodic
   *  \warning untested
   */ 
  static ICoord getNeighborIndex(ICoord index, int level, Direction dir)
  {   
    ICoord iLen   = voxelIntLength(level);
    int dec       = 0;    
    ICoord index2 = 0;
    ICoord coord = (index&INDEX_MASK) + getOrientation(dir,0)*iLen;

    if (BT != BoundaryType::PERIODIC)
      {
	if ((coord&(BBOX_ILEN-1)) != coord)
	  return index;
	index2 = coord;
      }
    else index2 = (coord&(BBOX_ILEN-1)); // periodic boundary  
    
    for (int i=1;i<NDIM;++i)
      {
	index>>= INDEX_DEC;
	dec   += INDEX_DEC;
	ICoord coord = (index&INDEX_MASK) + getOrientation(dir,i)*iLen;
	
	if (BT != BoundaryType::PERIODIC)
	  {
	    if ((coord&(BBOX_ILEN-1)) != coord)
	      return index;
	    index2 |= (coord<<dec);
	  }
	else index2 |= (coord&(BBOX_ILEN-1))<<dec; // periodic boundary
	  	  
      }
    return index2;
  }

   /** \brief Get the index of the same level direct neighbor of a voxel with index 'index'
   *  and level 'level'. Note that these voxel do not have to exist. The boundary conditions
   *  are considered to be periodic whatever the configuration of the grid.
   *  \param index the index whose neighbor we look for 
   *  \param level the level of the current index
   *  \param dir the direction toward the neighbor
   *  \return the index of the neighbor 
   *  \warning untested
   */ 
  static ICoord getNeighborIndex_periodic(ICoord index, int level, Direction dir)
  {   
    ICoord iLen   = voxelIntLength(level);
    int dec       = 0;    
    ICoord coord = (index&INDEX_MASK) + getOrientation(dir,0)*iLen;
    ICoord index2 = (coord&(BBOX_ILEN-1)); // periodic boundary  
    
    for (int i=1;i<NDIM;++i)
      {
	index>>= INDEX_DEC;
	dec   += INDEX_DEC;
	ICoord coord = (index&INDEX_MASK) + getOrientation(dir,i)*iLen;
	index2 |= (coord&(BBOX_ILEN-1))<<dec; // periodic boundary	  	  
      }
    return index2;
  }

   /*  
  //! gets the index of a voxel at level MAX_LEVEL_FROM_ROOT containing the point with 
  //! coordinates 'coords'.
  // FIXME : use geometry for boundary conditions !!!!!
  template <class InputIterator>
  ICoord coords2Index(InputIterator coords)
  {
    ICoord index=0; 
    ICoord c;

    for (int i=0;i<NDIM;++i)
      {
	// c is the integer coord of the center of the voxel @MAX_LEVEL
	double dc=geometry->checkBoundary(*coords,i);
	if (dc<=xMin[i]) c=1;
	else if (dc>=xMax[i]) c=BBOX_ILEN-1;
	else c=ICoord(((dc - xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);	
	index += (c<<(INDEX_DEC*i));
	++coords;
      }

    return index;
  }
*/
  
  /** \brief Compute the 'index' of a voxel at level 'level' (or MAX_MAX_LEVEL_FROM_ROOT by default)
   *  containing the 
   *  point with coordinates coords. Note that the voxel does NOT have to exist in
   *  the AMR grid. Moreover, this will always return a valid index, even when the coords
   *  are out of boundary and boundary conditions are non periodic !
   * \param coords an iterator to the coordinates of the point
   * \param level compute the index a voxel at level 'level' would have. By default, the index
   *  at maximum level is computed.
   * \tparam CT the type of 
   * \return the index of the voxel.
   * \see coords2index_noCheck
   */ 
  template <class InputIterator>
  ICoord coords2Index(InputIterator coords, int level=MAX_LEVEL_FROM_ROOT) const
  {
    static const ICoord BBOX_ILEN_CPY=BBOX_ILEN;
    ICoord c; // c is the integer coord of the center of the voxel @MAX_LEVEL
    ICoord dec=0;
    auto dc=geometry->checkBoundary(*coords,0);
    //double dc=geometry->checkBoundary(*coords,0);

    if (dc<xMin[0]) c=1;
    else if (dc>=xMax[0]-halfEpsilon[0]) c=BBOX_ILEN-1;
    //else c=ICoord(((dc - xMin[0])*(deltaX_inv[0])) * BBOX_ILEN) | ICoord(1);
    else c = hlp::numericStaticCast<ICoord>
	   (((dc - xMin[0])*(deltaX_inv[0])) * BBOX_ILEN_CPY) | ICoord(1);
	   
    //else c=ICoord(((dc+halfEpsilon[0] - xMin[0])*(deltaX_inv[0])) * BBOX_ILEN) | ICoord(1);

    ICoord index=c; 
    for (int i=1;i<NDIM;++i)
      {
	dec+=INDEX_DEC;
	++coords;

	dc=geometry->checkBoundary(*coords,i);
	if (dc<xMin[i]) c=1;
	else if (dc>=xMax[i]-halfEpsilon[i]) c=BBOX_ILEN-1;
	//else c=ICoord(((dc - xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);
	else c = hlp::numericStaticCast<ICoord>
	       (((dc - xMin[i])*(deltaX_inv[i])) * BBOX_ILEN_CPY) | ICoord(1);
	//else c=ICoord(((dc+halfEpsilon[i] - xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);	
	index |= (c<<dec);	
      }

    if (level!=MAX_LEVEL_FROM_ROOT) index = setIndexLevel(index,level);

    return index;
  }

  /** \brief Compute the 'index' of a voxel at MAX_LEVEL_FROM_ROOT containing the 
   *  point with coordinates coords. Note that the voxel does NOT have to exist in
   *  the AMR grid. Moreover, this will always return a valid index, even when the coords
   *  are out of boundary and boundary conditions are non periodic !
   * \param coords an iterator to the coordinates of the point
   * \param warning a flag whose nth bit is set to 1 if the nth coord is out of bound.
   * \param includeRightBoundary if true, then issue a wrning also if the point is exactly
   * on the right boundary (i.e. according to simulation of simplicity convention, the 
   * right boundary should NOT be within the box for non periodic boundary conditions)
   * \return the index of the voxel.
   * \see coords2index , coords2index_noCheck
   */ 
  template <class InputIterator, class T>
  ICoord coords2Index(InputIterator coords, T &warning, bool includeRightBoundary=true) const
  {
    static const ICoord BBOX_ILEN_CPY=BBOX_ILEN;
    ICoord c; // c is the integer coord of the center of the voxel @MAX_LEVEL
    ICoord dec=0;
    auto dc=geometry->checkBoundary(*coords,0);

    
    //double dc=geometry->checkBoundary(*coords,0);
    warning=0;
    if (dc<xMin[0]) {c=1;warning|=(1L<<0);}
    else if (dc>=xMax[0]-halfEpsilon[0]) 
      {
	c=BBOX_ILEN-1;
	if ((!includeRightBoundary)||(dc>xMax[0])) 
	  warning|=(1L<<0);
      }
    else 
      {
	auto val = (dc-xMin[0])*(deltaX_inv[0])*BBOX_ILEN_CPY;
	//double val = (dc-xMin[0])*(deltaX_inv[0])*BBOX_ILEN;
#ifdef HAVE_BOOST
#ifdef HAVE_GMP
	if ( fabs(val-trunc(val)) < halfEpsilonILen[0])
	  {
	    //typedef boost::multiprecision::mpf_float mpfloat;
	    
	    mpfloat mp_val = 
	      hlp::numericStaticCast<mpfloat>(geometry->checkBoundary(*coords,0));
	    //mpfloat mp_val = hlp::numericStaticCast<mpfloat>(dc);
	    mp_val-=mp_xMin[0];
	    mp_val*=mp_BBOX_ILEN;
	    mp_val/=mp_deltaX[0];
	    c = mp_val.convert_to<ICoord>() | ICoord(1) ;	    
	  }
	else 
#endif
#endif
	  c=hlp::numericStaticCast<ICoord>(val)|ICoord(1);
	//c=static_cast<ICoord>(val)|ICoord(1);
	//c=ICoord(val) | ICoord(1);
      }
    //c=ICoord(((dc+halfEpsilon[0]- xMin[0])*(deltaX_inv[0])) * BBOX_ILEN) | ICoord(1);

    ICoord index=c; 
    for (int i=1;i<NDIM;++i)
      {
	dec+=INDEX_DEC;
	++coords;	
	dc=geometry->checkBoundary(*coords,i);
	if (dc<xMin[i]) {c=1;warning|=(1L<<i);}
	else if (dc>=xMax[i]-halfEpsilon[i]) 
	  {
	    c=BBOX_ILEN-1;
	    if (!includeRightBoundary||(dc>xMax[i])) warning|=(1L<<i);
	  }
	else 
	  {
	    auto val = (dc-xMin[i])*(deltaX_inv[i])*BBOX_ILEN_CPY;
#ifdef HAVE_BOOST
#ifdef HAVE_GMP
	    if ( fabs(val-trunc(val)) < halfEpsilonILen[i])
	      {
		//typedef boost::multiprecision::mpf_float mpfloat;
		mpfloat mp_val = 
		  hlp::numericStaticCast<mpfloat>(geometry->checkBoundary(*coords,i));
		//mpfloat mp_val = hlp::numericStaticCast<mpfloat>(dc);
		mp_val-=mp_xMin[i];
		mp_val*=mp_BBOX_ILEN;
		mp_val/=mp_deltaX[i];
		c = mp_val.convert_to<ICoord>() | ICoord(1);
	      }
	    else 
#endif
#endif
	      c=hlp::numericStaticCast<ICoord>(val)|ICoord(1);
	    //c=static_cast<ICoord>(val)|ICoord(1);
	    //c=ICoord(val) | ICoord(1);
	  }
	//c=ICoord(((dc +halfEpsilon[i]- xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);
	index |= (c<<dec);	
      }
    
    return index;
    
    /*
    ICoord c; // c is the integer coord of the center of the voxel @MAX_LEVEL
    ICoord dec=0;
    double dc=geometry->checkBoundary(*coords,0);
    warning=0;
    if (dc<xMin[0]) {c=1;warning|=(1L<<0);}
    else if (dc>=xMax[0]-halfEpsilon[0]) 
      {
	c=BBOX_ILEN-1;
	if ((!includeRightBoundary)||(dc>xMax[0])) 
	  warning|=(1L<<0);
      }
    else c=ICoord(((dc+halfEpsilon[0]- xMin[0])*(deltaX_inv[0])) * BBOX_ILEN) | ICoord(1);

    ICoord index=c; 
    for (int i=1;i<NDIM;++i)
      {
	dec+=INDEX_DEC;
	++coords;

	dc=geometry->checkBoundary(*coords,i);
	if (dc<xMin[i]) {c=1;warning|=(1L<<i);}
	else if (dc>=xMax[i]-halfEpsilon[i]) 
	  {
	    c=BBOX_ILEN-1;
	    if (!includeRightBoundary||(dc>xMax[i])) warning|=(1L<<i);
	  }
	else c=ICoord(((dc +halfEpsilon[i]- xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);	
	index |= (c<<dec);	
      }
    
    return index;
    */
  }

  /** \brief Compute the 'index' of a voxel at MAX_LEVEL_FROM_ROOT containing the 
   *  point with coordinates coords without performing any test on the validity of coords.
   *  Note that the voxel does NOT have to exist in
   *  the AMR grid.
   * \param coords an iterator to the coordinates of the point
   * \return the index of the voxel.
   * \see coords2index
   */ 
  
  template <class InputIterator>
  ICoord coords2Index_noCheck(InputIterator coords)
  {
    return coords2Index(coords);
    ICoord c=ICoord(((coords[0] - xMin[0])*(deltaX_inv[0])) * BBOX_ILEN) | ICoord(1);
    ICoord dec=0;

    ICoord index=c; 
    for (int i=1;i<NDIM;++i)
      {
	dec+=INDEX_DEC;
	++coords;
	c=ICoord(((coords[i] - xMin[i])*(deltaX_inv[i])) * BBOX_ILEN) | ICoord(1);
	index |= (c<<dec);	
      }

    return index;
  }
  
  /*
  //! gets the coordinates of the center of a voxel with index index
  template <class OutputIterator>
  void index2Coords(ICoord index, OutputIterator coords)
  {
    for (int i=0;i<NDIM;++i)
      {
	(*coords)=xMin[i]+ICoordUnitLen[i]*((index>>(i*INDEX_DEC))&INDEX_MASK);
	++coords;
      }
  }
  */

  /** \brief Compute the integer coordinates of the center of a voxel given its index.
   * \param index the index of the voxel
   * \param level the level of the retrieved iCoords.
   * \param[out] iCoords the integer coordinates along each dimension   
   * \warning untested
   */
  static void index2ICoords(ICoord index, ICoord iCoords[NDIM],
			    int level=MAX_LEVEL_FROM_ROOT)
  {
    const int levelDec = (MAX_LEVEL_FROM_ROOT - level + 1);
    //const ICoord hl = voxelIntHalfLength(level);

    for (int i=0;i<NDIM;++i)
      {
	//iCoords[i]=((index&(INDEX_MASK))>>levelDec) + hl;
	iCoords[i]=((index&INDEX_MASK)>>levelDec);
	index >>= INDEX_DEC;
      }
  }
  

   /** \brief Compute the index of a voxel at a given level and with given integer 
    *  coordinates
    * \param iCoords the integer coordinates of the voxel lower left corner o
    * \param level the level of the voxel 
    * \warning untested
    */
  /*
  static ICoord iCoords2Index(ICoord iCoords[NDIM], int level)
  {
    ICoord index=0;
    int dec=0;
    const int levelDec = (MAX_LEVEL_FROM_ROOT - level + 1);
    const ICoord hl = voxelIntHalfLength(level);
    
    for (int i=0;i<NDIM;++i)
      {
	index|=((iCoords[i]>>levelDec)+hl)<<dec;
	dec+=INDEX_DEC;
      }
    return index;
  }
  */

  /** \brief Compute the coordinates of the center of a voxel with index 'index'
   * \param index the index of the voxel
   * \param[out] coords an output iterator to the coordinates of the point
   */
  template <class OutputIterator>
  void index2CenterCoords(ICoord index, OutputIterator coords) const 
  {
    (*coords)=xMin[0]+ICoordUnitLen[0]*(index&INDEX_MASK);
    for (int i=1;i<NDIM;++i)
      {
	index>>=INDEX_DEC;
	++coords;
	(*coords)=xMin[i]+ICoordUnitLen[i]*(index&INDEX_MASK);	
      }
  }

  // FIXME : Could be faster ????
  /** \brief compute the coordinates of a corner vertex of a voxel with index 'index' at
   *  level 'level'.
   * Which vertex is indicated by the bits sets in corner. 
   * e.g. in 3D :
   * corner = 000b indicates left most corner in each dimension 
   * corner = 011b indicates vertex with coordinates (x,y,z)=[1,1,0]
   * corner = 111b indicates the right most corner in each dimension
   * \param index the index of the voxel
   * \param level the level of the voxel
   * \param[out] coords an output iterator to the coordinates of the corner vertex
   * \param corner An index representing which corner we want to compute
   */  
  template <class OutputIterator>
  void index2CornerCoords(ICoord index, int level, OutputIterator coords, int corner=0) const
  {
    //double hl = getVoxelHalfLength(level,0)*(((corner&1)<<1)-1);
    long hl = long(voxelIntHalfLength(level))*(((corner&1)<<1)-1);
    (*coords)=xMin[0]+ICoordUnitLen[0]*( (index&INDEX_MASK) + hl);//+hl;  
    for (int i=1;i<NDIM;++i)
      {
	index>>=INDEX_DEC;
	corner>>=1;
	//hl = getVoxelHalfLength(level,i)*(((corner&1)<<1)-1);
	hl = long(voxelIntHalfLength(level))*(((corner&1)<<1)-1);
	++coords;
	(*coords)=xMin[i]+ICoordUnitLen[i]*( (index&INDEX_MASK) + hl);//+hl;	
      }
  }

  /** \brief computes the coordinates of a corner vertex of a voxel with index 'index' at
   * level 'level' as well as the coordinates of its opposite vertex in the voxel.
   * Which vertex is indicated by the bits sets in corner. This version also retrieves
   * the coordinates of the opposite vertex (the two give a bounding box)
   * e.g. in 3D :
   * corner = 000b indicates left most corner in each dimension 
   * corner = 011b indicates vertex with coordinates (x,y,z)=[1,1,0]
   * corner = 111b indicates the right most corner in each dimension
   * \param index the index of the voxel
   * \param level the level of the voxel
   * \param[out] coords an output iterator to the coordinates of the corner vertex
   * \param[out] oppCoords an output iterator to the coordinates of the opposite
   * corner vertex
   * \param corner An index representing which corner we want to compute
   */
  template <class OutputIterator>
  void index2CornerCoordsAndOpp(ICoord index, int level, 
				OutputIterator coords, 
				OutputIterator oppCoords,
				int corner=0) const
  {
    //ICoord levelMask = (~((1L<<(MAX_LEVEL_FROM_ROOT - level))-1)) & INDEX_MASK;

    long hl = long(voxelIntHalfLength(level))*(((corner&1)<<1)-1);

    (*coords)=xMin[0]+ICoordUnitLen[0]*( (index&INDEX_MASK) + hl);
    (*oppCoords)=xMin[0]+ICoordUnitLen[0]*( (index&INDEX_MASK) - hl);  
    
    for (int i=1;i<NDIM;++i)
      {
	index>>=INDEX_DEC;
	corner>>=1;
	hl = long(voxelIntHalfLength(level))*(((corner&1)<<1)-1);
	++coords;
	++oppCoords;

	(*coords)=xMin[i] + ICoordUnitLen[i]*( (index&INDEX_MASK) + hl );
	(*oppCoords)=xMin[i] + ICoordUnitLen[i]*( (index&INDEX_MASK) - hl ); 
      }

    /*
    double hl = getVoxelHalfLength(level,0)*(((corner&1)<<1)-1);    
    double center = xMin[0]+ICoordUnitLen[0]*(index&INDEX_MASK);
    (*coords)=center+hl;
    (*oppCoords)=center-hl;
    for (int i=1;i<NDIM;++i)
      {
	index>>=INDEX_DEC;
	corner>>=1;
	hl = getVoxelHalfLength(level,i)*(((corner&1)<<1)-1);
	++coords;
	++oppCoords;
	center = xMin[i]+ICoordUnitLen[i]*(index&INDEX_MASK);
	(*coords)=center+hl;
	(*oppCoords)=center-hl;
      }
    */
  }

  /** \brief given an index at any level, returns the index of a voxel 
   *  at level 'level' that would contain the center of the voxel it represents.
   *  Note that 'index' does not have to represent a valid voxel in the grid.
   *  \return the new index at level 'level'
   */
  static ICoord setIndexLevel(ICoord index, int level)
  {   
    const int levelDec = (MAX_LEVEL_FROM_ROOT - level + 1);
    const ICoord hl = voxelIntHalfLength(level);
    ICoord factor = (index & INDEX_MASK)>>levelDec;
    ICoord result = ((factor<<levelDec) + hl);
    int dec = INDEX_DEC;
    for (int i=1;i<NDIM;++i)
      {	
	index  >>= INDEX_DEC;
	factor = (index & INDEX_MASK)>>levelDec;
	result |= (((factor<<levelDec) + hl)<<dec) ;
	dec+=INDEX_DEC;
      }
    return result;
  }

  /** \brief check if a voxel is on the boundary of the bounding box along each dimension
   * \param index the index of the voxel to check
   * \param level the level of the voxel
   * \param onBoundary an array where the result is stored for each dimension. The value is
   *  set to {-1,+1,0} if the voxel is on the left boundary, on the right boundary or not
   *  on the boundary along each dimension respectively.
   * \return true if the voxel is on any boundary, false otherwise: the ith bit of the 
   *  result is set if the voxel is on the corresponding dimension boundary.
   */
  template <class T>
  int isOnBoundary(ICoord index, int level, T onBoundary[NDIM])
  {
    int result=0;
    ICoord hlen=voxelIntHalfLength(level);
    
    for (int i=0;i<NDIM;++i)
      {
	ICoord coord= index&INDEX_MASK;
    
	if (coord<=hlen) {onBoundary[i]=-1;result&=(1<<i);}
	else if (coord>=BBOX_ILEN-hlen) {onBoundary[i]=1;result&=(1<<i);}
	coord>>=INDEX_DEC;
      }

    return result;
  }

  /*
  //! deduce the index of the child with rank "quadrant" from the index of the current
  //! voxel 'index' and its level 'level'. We could deduce 'level' from 'index', but that
  //! would waste computations ...
  static ICoord computeChildIndex(ICoord index, int level, int quadrant)
  {
    ICoord hl = voxelIntHalfLength(level+1);
    ICoord result=index;
    for (int i=0;i<NDIM;++i)
      {	
	if (quadrant&(1L<<i))
	  result += hl;
	else
	  result -= hl;
	hl<<=INDEX_DEC;
      }
    return result;
  } 
  */

  
  /** \brief deduce the index of the child in the quadrant "quadrant" of the voxel with
   * index 'index' and at level 'level'.
   * We could deduce 'level' from 'index', but that would waste computations ... 
   * \param index the index of the voxel
   * \param level the level of the voxel
   * \param quadrant the quadrant of the voxel the child falls into
   * \return the index of the child voxel.
   */  
  static ICoord computeChildIndex(ICoord index, int level, int quadrant)
  {
    ICoord hl = voxelIntHalfLength(level+1);
    ICoord result= index + hl*(((quadrant&1)<<1)-1);
    for (int i=1;i<NDIM;++i)
      {
	quadrant>>=1;
	hl<<=INDEX_DEC;
	result += hl*(((quadrant&1)<<1)-1);	
      }
   
    return result;
  } 
  
  /* 
  // FIXME: This version may be safer, check it (-> possible overflow for signed ICoord)
  static ICoord computeChildIndex(ICoord index, int level, int quadrant)
  {
    ICoord hl = voxelIntHalfLength(level+1);
    ICoord result= (index&INDEX_MASK) + hl*(((quadrant&1)<<1)-1);
    ICoord dec = INDEX_DEC;
    for (int i=1;i<NDIM;++i)
      {
	quadrant>>=1;
	index>>=INDEX_DEC;
	result |= ((index&INDEX_MASK) + hl*(((quadrant&1)<<1)-1))<<dec;
	dec += INDEX_DEC;
      }
   
    return result;
  } 
  */
  
  

  /** \brief retrieve the (1<<NDIM) leaves sharing the vertex \a vertId of voxel. Note that 
   *  the first returned voxel is always \a voxel even when it is not a leaf. Note also that
   *  a voxel at lower level than \a voxel may be replicated in the output when the vertex
   *  falls on one of its faces.
   * \param voxel A pointer to the reference voxel
   * \param vertexId Index of the vertex
   * \param[out] out An output iterator to Voxel pointers where the (1<<NDIM) neighbor 
   * voxels will be stored. In case a given neighbor does not exist (non periodic boundary 
   * conditions), a pointer to the queried \a voxel is stored instead.
   * \return the number of non-existing neighbors (due to non periodic boundary conditions)
   */ 
  template <class OutputIterator>
  int getVertexNeighbors(Voxel *voxel, int vertexId, OutputIterator out)
  {
    long iLen = voxelIntHalfLength(voxel->getLevel())-1;
    ICoord index = voxel->getIndex();    
    //ICoord iDelta[NDIM];
    int count = 0; // number of neighbors that do not exist !

    (*out)=voxel;++out;
    for (int j=1;j<(1L<<NDIM);++j)
      {
	ICoord id=index;	
	ICoord vid=vertexId;
	ICoord newIndex=0;
	ICoord dec=0;
	Voxel *result=NULL;

	for (int i=0;i<NDIM;++i)
	  {	
	    // ((j>>i)&1) -> should we move along dimension i ?
	    // (((vid&1)<<1)-1) ->  orientation of the movement (+1/-1)
	    
	    //ICoord coord = (id&INDEX_MASK) + ((j>>i)&1)*(((vid&1)<<1)-1)*iLen;
	    if (BT != BoundaryType::PERIODIC)
	      {
		long mv= ((j>>i)&1)<<1; // should we move along dimension i  (0/2)?
		long dir= ((long(vid&1)<<1)-1); // orientation of the movement (+1/-1)
		long coord = long(id&INDEX_MASK) + dir*(iLen + mv);
	    
		if ((coord&(BBOX_ILEN-1)) != coord) 
		  result=voxel; //This one is out of bound !
		else newIndex |= (ICoord(coord)<<dec);
	      }
	    else //newIndex |= (coord&(BBOX_ILEN-1))<<dec; // periodic boundary
	      {
		long mv= ((j>>i)&1)<<1; // should we move along dimension i  (0/2)?
		long dir= ((long(vid&1)<<1)-1); // orientation of the movement (+1/-1)
		long coord = long( (id&INDEX_MASK) + BBOX_ILEN) + dir*(iLen + mv);
	
		newIndex |= ((ICoord(coord)&(BBOX_ILEN-1))<<dec);
	      }
	    dec+=INDEX_DEC;
	    id>>=INDEX_DEC;
	    vid>>=1;
	  }
	
	if (result!=voxel) 
	  result=getVoxel(newIndex);
	else count++;
	
	(*out)=result;
	++out;
      }

    return count;
  }

  /** \brief get the direction along each axis toward each neighbor of a vertex. This is
   *  usefull to retrieve the position of getVertexNeighbors() output.  
   * \param vertexId Index of the vertex
   * \param[out] direction The direction (+1/-1) along each of the NDIM dimension of each 
   * of the (1<<NDIM) neighbors.
   */ 
  template <class T>
  void getVertexNeighborsDirection(int vertexId, T (&direction)[1L<<NDIM][NDIM])				   
  {
    for (int j=0;j<(1L<<NDIM);++j)
      {
	int vid=vertexId;
	for (int i=0;i<NDIM;++i)
	  {
	    // ((j>>i)&1) -> true if we should move along dimension i
	    // (vid&1) ->  true if we should move in the positive direction
	    // ((((j>>i)&1)<<1)-1) = {+1,-1} if we {should,should not} move along i
	    // (((vid&1)<<1)-1) = {+1,-1} if we move in the {positive,negative} direction
	    direction[j][i]=((((j>>i)&1)<<1)-1)*(((vid&1)<<1)-1);
	    vid>>=1;
	  }
      }
  }

   /** \brief get the direction along each axis toward the neighbor with index 'index' of 
    * a vertex. This is
    *  usefull to retrieve the position of getVertexNeighbors() output.  
    * \param vertexId Index of the vertex
    * \param index the index of the neighbor
    * \param[out] direction The direction (+1/-1) along each of the NDIM dimension of each 
    * of the (1<<NDIM) neighbors.
    */ 
  template <class T>
  void getVertexNeighborDirection(int vertexId, int index, T (&direction)[NDIM])				   
  {    

    for (int i=0;i<NDIM;++i)
      {
	// ((index>>i)&1) -> true if we should move along dimension i
	// (vid&1) ->  true if we should move in the positive direction
	// ((((index>>i)&1)<<1)-1) = {+1,-1} if we {should,should not} move along i
	// (((vid&1)<<1)-1) = {+1,-1} if we move in the {positive,negative} direction
	direction[i]=((((index>>i)&1)<<1)-1)*(((vertexId&1)<<1)-1);
	vertexId>>=1;
      }     
  }

  /** \brief retrieve the (1<<(NDIM-1)) leaves sharing the segment adjacent to vertex 
   *  \a vertId and along direction dir. The first returned voxel is always \a voxel 
   *  even when it is not a leaf. Note also that a voxel at lower level than \a voxel 
   *  may be replicated in the output when the vertex falls on one of its faces.
   * \param voxel A pointer to the reference voxel
   * \param vertexId Index of the vertex
   * \param dim dimension along which the segment is oriented
   * \param[out] out An output iterator to Voxel pointers where the (1<<NDIM) neighbor 
   * voxels will be stored. In case a given neighbor does not exist (non periodic boundary 
   * conditions), a pointer to the queried \a voxel is stored instead.
   * \return the number of non-existing neighbors (due to non periodic boundary conditions)
   */ 
  template <class OutputIterator>
  int getSegmentNeighbors(Voxel *voxel, int vertexId, int dim, OutputIterator out)
  {
    long ddim = 1L<<dim;
    long iLen = voxelIntHalfLength(voxel->getLevel())-1;
    ICoord index = voxel->getIndex();    
    //ICoord iDelta[NDIM];
    int count = 0; // number of neighbors that do not exist !

    (*out)=voxel;++out;
    for (int j=1;j<(1L<<NDIM);++j)
      {
	if (j&ddim) continue;//j+=ddim;

	ICoord id=index;	
	ICoord vid=vertexId;
	ICoord newIndex=0;
	ICoord dec=0;
	Voxel *result=NULL;

	for (int i=0;i<NDIM;++i)
	  {	
	    // ((j>>i)&1) -> should we move along dimension i ?
	    // (((vid&1)<<1)-1) ->  orientation of the movement (+1/-1)
	    
	    //ICoord coord = (id&INDEX_MASK) + ((j>>i)&1)*(((vid&1)<<1)-1)*iLen;
	    if (BT != BoundaryType::PERIODIC)
	      {
		long mv= ((j>>i)&1)<<1; // should we move along dimension i  (0/2)?
		long dir= ((long(vid&1)<<1)-1); // orientation of the movement (+1/-1)
		long coord = long(id&INDEX_MASK) + dir*(iLen + mv);
	    
		if ((coord&(BBOX_ILEN-1)) != coord) 
		  result=voxel; //This one is out of bound !
		else newIndex |= (ICoord(coord)<<dec);
	      }
	    else //newIndex |= (coord&(BBOX_ILEN-1))<<dec; // periodic boundary
	      {
		long mv= ((j>>i)&1)<<1; // should we move along dimension i  (0/2)?
		long dir= ((long(vid&1)<<1)-1); // orientation of the movement (+1/-1)
		long coord = long( (id&INDEX_MASK) + BBOX_ILEN) + dir*(iLen + mv);
	
		newIndex |= ((ICoord(coord)&(BBOX_ILEN-1))<<dec);
	      }
	    dec+=INDEX_DEC;
	    id>>=INDEX_DEC;
	    vid>>=1;
	  }
	
	if (result!=voxel) 
	  result=getVoxel(newIndex);
	else count++;
	
	(*out)=result;
	++out;
      }

    return count;
  }

  /** \brief get the direction along each axis toward each neighbor of a segment. This is
   *  usefull to retrieve the position of getSegmentNeighbors() output.  
   * \param vertexId Index of the vertex
   * \param dim dimension along which the segment is oriented
   * \param[out] direction The direction (+1/-1) along each of the NDIM dimension of each 
   * of the (1<<(NDIM-1)) neighbors.
   */ 
  template <class T>
  void getSegmentNeighborsDirection(int vertexId, int dim, T direction[1L<<(NDIM-1)][NDIM])				   
  {
    long ddim = 1L<<dim;
    int k=0;
    for (int j=0;j<(1L<<NDIM);++j)
      {
	if (j&ddim) continue;//j+=ddim;

	int vid=vertexId;
	for (int i=0;i<NDIM;++i)
	  {	    
	    // ((j>>i)&1) -> true if we should move along dimension i
	    // (vid&1) ->  true if we should move in the positive direction
	    // ((((j>>i)&1)<<1)-1) = {+1,-1} if we {should,should not} move along i
	    // (((vid&1)<<1)-1) = {+1,-1} if we move in the {positive,negative} direction
	    direction[k][i]=((((j>>i)&1)<<1)-1)*(((vid&1)<<1)-1);
	    vid>>=1;
	  }
	k++;
	// By convention, direction is positive along the segment
	// direction[k++][dim]=1;
      }
  }

  /** \brief Retrieve all the leaf voxels that would overlap a box defined by its 
   *  coordinates and dimensions within a uniform grid of given resolution. The resolution
   *  of the uniform grid is defined by that of a voxel at level \a 'level'. Note that
   *  level may take negative values (i.e. when the uniform grid has lower resolution than
   *  the root voxels of the AMR grid).   
   *  \param position The coordinates of the lower left corner of the box within 
   *   the uniform grid.
   *  \param dims The dimensions of the box in units of the uniform grid voxels.
   *  \param level The resolution of the uniform grid. This may be negative.
   *  \param[out] out An output iterator where the intersecting leaf voxel will be stored
   *  \param nonNullOnly if true, only retrieve the leaves with non null data.
   *  \tparam T The coordinate type of the box
   *  \tparam OutputIterator The type of the output iterator \a out
   */
  template <class T, class OutputIterator>
  void getLeavesGridOverlap(const T position[NDIM], 
			    const T dims[NDIM], 
			    int level,
			    OutputIterator out,
			    bool nonNullOnly=false)
  {    
    internal::localAmrGridOverlapVisitor::
      GetLeavesBoxOverlapT<MyType,OutputIterator> visitor(out,nonNullOnly);

    visitGridOverlap(position,dims,level,visitor,nonNullOnly);
  }
  

  /** \brief Recursively visit the voxels that overlap a box defined by its 
   *  coordinates and dimensions within a uniform grid of given resolution. The resolution
   *  of the uniform grid is defined by that of a voxel at level \a 'level'. Note that
   *  level may take negative values (i.e. when the uniform grid has lower resolution than
   *  the root voxels of the AMR grid). The Visitor should implement a visit(Voxel *v) 
   *  and visited(Voxel *v) method that will be called in a way equivalent to :
   *
   *  \verbatim
   *  template <class V>
   *  void visitTree(Voxel *voxel, V &v)
   *   {
   *     if (v.visit(voxel))
   *     {
   *	  for (int i=0;i<nChildren;++i)
   *	    {  
   *          if (box.overlaps(voxel->getChild(i))) 
   *	        visitTree(voxel->getChild(i),v);
   *	    }
   *	  v.visited(voxel);
   *    }    
   *  }
   *
   *  foreach (rootVoxel *voxel)
   *    visitTree(voxel, visitor);
   *  \endverbatim
   *
   * where voxel is always a voxel that intersects or is inside the box. 
   *
   *
   *  \param position The integral coordinates of the lower left corner of the box within 
   *   the uniform grid in units of the root grid voxel
   *  \param dims The integral dimensions of the box in units of the root grid voxels.
   *  \param level The resolution of the uniform grid. This may be negative.
   *  \param v  the visitor (see description)
   *  \param nonNullOnly if true, only retrieve the leaves with non null data.
   *  \tparam T The coordinate type of the box
   *  \tparam Visitor The visitor type implementing a visit(Voxel *v) and a 
   *   visited(Voxel *v) method.
   *  \warning The bounding box defined by position and dims may not overlap the boundary
   */
  template <class T, class Visitor>
  void visitGridOverlap(const T position[NDIM], 
			const T dims[NDIM], 
			int level,
			Visitor &v,
			bool nonNullOnly=false)
  {
    internal::localAmrGridVisitor::BoxOverlapHandlerT<MyType,Visitor> visitor(v);
    ICoord rootMin[NDIM];
    ICoord rootMax[NDIM];
    ICoord rootDims[NDIM];

    ICoord iBox[2][NDIM];
    ICoord voxIDim[NDIM+1]; 
    ICoord iStride[NDIM];

    std::copy(position,position+NDIM,rootMin);
    std::copy(dims,dims+NDIM,rootDims);
    for (int i=0;i<NDIM;++i)
      rootMax[i]=rootMin[i]+rootDims[i]-1;

    if (level<0)
      {
	for (int j=0;j<NDIM;++j)
	  {
	    iBox[0][j]=rootMin[j]*(ROOT_VOXEL_ILEN<<(-level));
	    iBox[1][j]=rootMax[j]*(ROOT_VOXEL_ILEN<<(-level));
	    rootMin[j]<<=(-level);
	    rootDims[j]<<=(-level);
	    rootMax[j]=rootMin[j]+rootDims[j]-1;	    
	  }
      }
    else if (level>0)
      {
	for (int j=0;j<NDIM;++j)
	  {
	    iBox[0][j]=rootMin[j]*(ROOT_VOXEL_ILEN>>(level));
	    iBox[1][j]=rootMax[j]*(ROOT_VOXEL_ILEN>>(level));
	    rootMin[j]>>=level;	      
	    rootMax[j]>>=level;
	    rootDims[j]=rootMax[j]-rootMin[j]+1;
	  }
      }          
    
    ICoord rootIndex=0;
    for (int i=0;i<NDIM;++i)
      {
	rootIndex += rootMin[i]*rootIStride[i];		
	iStride[i] = rootIStride[i+1]-(rootDims[i]*rootIStride[i]);	  
      }

    voxIDim[NDIM]=ROOT_VOXEL_IHALFLEN;
      
    ICoord w[NDIM]={0};
    Voxel *voxel=&rootVoxels[rootIndex];
    for (;w[NDIM-1]<rootDims[NDIM-1];++voxel)
      {
	if ((!nonNullOnly)||(!voxel->isLeaf())||(voxel->data!=0))
	  {
	    for (int i=0;i<NDIM;++i)
	      voxIDim[i] = ROOT_VOXEL_ILEN*(rootMin[i]+w[i]);

	    visitor.initialize_overlap_handler(iBox,voxIDim,voxIDim[NDIM]);
	    visitTree_rec(voxel,visitor);
	  }

	hlp::getNext<NDIM>(w,voxel,rootDims,iStride);		 
      }    
  }


  /** \brief Retrieve all the leaf voxels that overlap a given bounding box. This function
   *   also works for periodic boundary conditions. In that case, the bounding box may 
   *   overlap the AMR grid boundaries, but its extent along each dimension may not
   *   be larger than half the AMR box size.
   *  \param box The box to intersect, where \a box[0] represent the lefmost coordinates 
   *  and \a box[1] The righmost coordinates.
   *  \param[out] out An output iterator where the intersecting leaf voxel will be stored
   *  \tparam T The coordinate type of the box
   *  \tparam OutputIterator The type of the output iterator \a out
   */
  template <class T, class OutputIterator>
  void getLeavesBBoxOverlap(const T (&box)[2][NDIM], OutputIterator out) 
  //, bool prt=false)
  {
    internal::localAmrGridOverlapVisitor::
      GetLeavesBoxOverlapT<MyType,OutputIterator> visitor(out);
    visitBBoxOverlap(box,visitor);
  }

  /** \brief Recursively visit the voxels that overlap a bounding box defined by minimal and
   *  maximal coordinates along each dimension. This function also works for periodic 
   *  boundary conditions and bounding boxes overlapping a boundary. 
   *  The Visitor should implement a visit(Voxel *v) and a 
   *  visited(Voxel *v) method that will be called in a way equivalent to :
   *
   *  \verbatim
   *  template <class V>
   *  void visitTree(Voxel *voxel, V &v)
   *   {
   *     if (v.visit(voxel))
   *     {
   *	  for (int i=0;i<nChildren;++i)
   *	    {  
   *          if (box.overlaps(voxel->getChild(i))) 
   *	        visitTree(voxel->getChild(i),v);
   *	    }
   *	  v.visited(voxel);
   *    }    
   *  }
   *
   *  foreach (rootVoxel *voxel)
   *    visitTree(voxel, visitor);
   *  \endverbatim   
   *
   * where voxel is always a voxel that intersects or is inside the box. 
   *   
   *  \param box The box to intersect, where \a box[0] represent the lefmost coordinates 
   *  and \a box[1] The righmost coordinates. 
   *  \param v  the visitor (see description)
   *  \tparam T The coordinate type of the box
   *  \tparam Visitor The visitor type implementing a visit(Voxel *v) and a 
   *   visited(Voxel *v) method.
   */
  template <class T, class Visitor>
  void visitBBoxOverlap(const T (&box)[2][NDIM], Visitor &v)
  {
    internal::localAmrGridVisitor::BoxOverlapHandlerT<MyType,Visitor> visitor(v);
    ICoord rootOrigin[NDIM]; // coordinates of the origin root cell intersecting box
    ICoord w[NDIM]; //coordinates in pixels in the root array
    ICoord dim[NDIM]; //size in terms of pixels of the root array
    ICoord rootIndex=0; // index in root array
    ICoord iStride[NDIM]; // used in peridoic case only, for hlp::getNext() 

    // The ICoord coordinates + half-Isize of the current root voxel
    ICoord voxIDim[NDIM+1]; 
    ICoord iBox[2][NDIM]; // The ICoord coordinates of the box
    // The ICoord coordinates of the intersection of iBox anf the boundary of the current root voxel
    //ICoord curIBox[2][NDIM]; 

    voxIDim[NDIM]=ROOT_VOXEL_IHALFLEN;    
        
    if (BT != BoundaryType::PERIODIC)
      {
	double pos;
	ICoord inside=1;
	ICoord i1;
	//iStride[0]=1;

	for (int i=0;i<NDIM;++i)
	  {	 
	    // This epsilon is because when we fall on a left side of a voxel, we want
	    // to include its left neighbor also ...
	    pos=(box[0][i]-xMin[i]) - halfEpsilon[i]; 
	    if (pos<0) pos=0;
	    else if (pos>=deltaX[i]) 
	      {
		// This epsilon is because we want to include the right BBox boundaries !
		if (pos>=deltaX[i]+halfEpsilon[i]) inside=0; 
		pos=deltaX[i]-halfEpsilon[i]; 		
	      }
	    rootOrigin[i]=N_ROOT_PER_DIM * (pos * deltaX_inv[i]);
	    iBox[0][i]=BBOX_ILEN * (pos * deltaX_inv[i]);
	    
	    // Out choice of SoS implies that when we fall on the right side exactly,
	    // we also include the right neighbor !
	    pos=(box[1][i]-xMin[i]) + halfEpsilon[i];//+ getEpsilon(i); 
	    if (pos<0) {pos=0;inside=0;}
	    else if (pos>=deltaX[i]) pos=deltaX[i]-halfEpsilon[i];
	    i1=N_ROOT_PER_DIM * (pos * deltaX_inv[i]);
	    iBox[1][i]=BBOX_ILEN * (pos * deltaX_inv[i]);	

	    if (iBox[0][i]>iBox[1][i])
	      {
		std::swap(iBox[0][i],iBox[1][i]);
		std::swap(rootOrigin[i],i1);
	      }

	    w[i]=0;
	    dim[i]=inside+(i1-rootOrigin[i]);
	    rootIndex += rootOrigin[i]*rootIStride[i];
	    iStride[i] = rootIStride[i+1]-(dim[i]*rootIStride[i]);	    
	  }

	if (inside)
	  {
	    Voxel *voxel=&rootVoxels[rootIndex];
	    for (;w[NDIM-1]<dim[NDIM-1];++voxel)
	      {
		for (int i=0;i<NDIM;++i)
		  voxIDim[i] = ROOT_VOXEL_ILEN*(rootOrigin[i]+w[i]);	
		/*
		 printf("%f<(%f,%f)<%f  // %f<(%f,%f)<%f \n",
		 	   xMin[0],box[0][0],box[1][0],xMax[0],
		 	   xMin[1],box[0][1],box[1][1],xMax[1]);
		 printf("w0=%ld+%ld, w1=%ld+%ld / dim = (%ld,%ld) -> index=%ld\n",w[0],rootOrigin[0],w[1],rootOrigin[1],dim[0],dim[1],rootIndex+std::distance(&rootVoxels[rootIndex],voxel));
		*/
		/*
		if (voxel->isLeaf())
		  {
		    if ((!nonNull)||(voxel->data!=0))
		      {
			(*out)=voxel;
			++out;
		      }		      
		  }
		else 
		*/
		
		visitor.initialize_overlap_handler(iBox,voxIDim,voxIDim[NDIM]);
		visitTree_rec(voxel,visitor);
		//getLeavesBoxOverlap_rec(iBox,voxIDim,voxel,out);

		// ++w[0];++voxel;
		// if (w[0]>=dim[0]) 
		hlp::getNext<NDIM>(w,voxel,dim,iStride);		 
	      }
	  }
      }
    else // periodic boundary case
      {
	double x0,x1,dx;
	ICoord i1;
	// if (prt) printf("BBOX: min=(%g,%g) max=(%g,%g)\n",
	// 		box[0][0],box[0][1],
	// 		box[1][0],box[1][1]);
	for (int i=0;i<NDIM;++i)
	  {
	    x0 = geometry->checkBoundary(box[0][i],i) - xMin[i];
	    x1 = geometry->checkBoundary(box[1][i],i) - xMin[i];
	    dx = geometry->correctCoordsDiff(x1-x0,i);	    
	    if (dx<0)  std::swap(x0,x1);
	    
	    // This epsilon is because when we fall on a left side of a voxel, we want
	    // to include its left neighbor also. This is important as we do not know
	    // which simplex the vertex would fall in ...
	    if (x0<halfEpsilon[i])
	      {		    
		x0 = deltaX[i]-halfEpsilon[i];// xMax[i]-xMin[i];		
	      }
	    else x0-=halfEpsilon[i];

	    if (x1>deltaX[i]-halfEpsilon[i])//+xMax[i]-xMin[i])
	      {		    
		x1 = halfEpsilon[i];		
	      }
	    else x1+=halfEpsilon[i];
	    
	    // x1+=getEpsilon(i); //??
	    // This epsilon is because when we fall on a left side of a voxel, we want
	    // to include its left neighbor also ...
	    

	    rootOrigin[i]=ICoord(N_ROOT_PER_DIM*(x0*deltaX_inv[i]));//&(N_ROOT_PER_DIM-1);
	    i1=ICoord(N_ROOT_PER_DIM*(x1*deltaX_inv[i]));//&(N_ROOT_PER_DIM-1);
	  
	    iBox[0][i]=ICoord(BBOX_ILEN*(x0*deltaX_inv[i]));//&(BBOX_ILEN-1);
	    iBox[1][i]=ICoord(BBOX_ILEN*(x1*deltaX_inv[i]));//&(BBOX_ILEN-1);
	    
	    w[i]=0;
	    if (i1<rootOrigin[i])	  
	      dim[i]=1+i1+(N_ROOT_PER_DIM-rootOrigin[i]);	     
	    else	   
	      dim[i]=1+(i1-rootOrigin[i]);

	    // if (prt) printf("rootOrigin[%d]=%ld (x0=%g, x1=%g)\n",i,rootOrigin[i],x0,x1);
	    
	    /*
	    dim[i]=(i1-rootOrigin[i]);	   
	    if (dim[i]>(N_ROOT_PER_DIM>>1)) 
	      dim[i]=1+(N_ROOT_PER_DIM-dim[i]);
	    else
	      dim[i]+=1;
	    */
	  }

	ICoord iBoxPer[2][NDIM];	
	for (;w[NDIM-1]<dim[NDIM-1];)
	  {	  
	    ICoord coord;
	    ICoord rootIndex=0;	    
	    // FIXME: need to check !!!!
	    for (int i=0;i<NDIM;++i)
	      {
		coord = ((rootOrigin[i]+w[i])&(N_ROOT_PER_DIM-1));
		rootIndex+=coord*rootIStride[i];
		voxIDim[i] = ROOT_VOXEL_ILEN*coord;

		if (iBox[1][i]>=iBox[0][i]) // no need to periodize
		  {
		    iBoxPer[0][i]=iBox[0][i];
		    iBoxPer[1][i]=iBox[1][i];		    
		  }
		else // we are crossing a boundary !
		  {
		    if (coord != (rootOrigin[i]+w[i])) // we are on iBox[1] side
		      {
			iBoxPer[0][i]=0;
			iBoxPer[1][i]=iBox[1][i];
		      }
		    else // we are on iBox[0] side
		      {
			iBoxPer[0][i]=iBox[0][i];
			iBoxPer[1][i]=(BBOX_ILEN-1);
		      }
		  }
	      }
	    // if (prt) 
	    //   {
	    // 	printf("x=[%ld %ld], iBoxPer=[%ld,%ld][%ld,%ld]\n",
	    // 	       w[0],w[1],
	    // 	       iBoxPer[0][0],iBoxPer[0][1],
	    // 	       iBoxPer[1][0],iBoxPer[1][1]);
	    //   }
	    Voxel *voxel=&rootVoxels[rootIndex];
	    // if (prt) voxel->print(this,"Testing ROOTVOXEL:");

	    /*
	    if (voxel->isLeaf())
	      {
		if ((!nonNull)||(voxel->data!=0))
		  {
		    (*out)=voxel;
		    ++out;
		  }		      
	      }
	    else
	      */

	    visitor.initialize_overlap_handler(iBoxPer,voxIDim,voxIDim[NDIM]);
	    visitTree_rec(voxel,visitor);
	    //getLeavesBoxOverlap_rec(iBoxPer,voxIDim,voxel,out);
	      
	    // ++w[0];
	    // if (w[0]>=dim[0]) 
	    hlp::getNext<NDIM>(w,dim);	    
	  }
      }

  }

  /** \brief Similar to visitTree(V & visitor, int nThreads) but attributes a different
   *  visitor to each thread to prevent locking when visitors private data need to 
   *  be updated.
   *
   *  Recursively visit the tree and call visitor.visit(Voxel *voxel) with the 
   *  current voxel as parameter when visiting a voxel for the first time, and 
   *  visitor.visited(Voxel *voxel) when all the descendants of a voxel have been visited.
   *  The descendants of a voxel \a voxel are visited only if visitor.visit(voxel) 
   *  returns true, and the branch is skipped otherwise. Note that if visitor.visit(voxel)
   *  returns false, then visitor.visited(voxel) will not be called. The visitor should 
   *  also implement a visitor.visit(Voxel *voxel,i) and visitor.visited(Voxel *voxel,i) 
   *  that act similarly, except that they are called before and after visiting the ith 
   *  child of a voxel respectively. The visitor will be used in a way equivalent to:
   *
   *  \verbatim
   *  template <class V>
   *  void visitTree_rec(Voxel *voxel, V &visitor)
   *   {
   *     if (visitor.visit(voxel))
   *     {
   *	  for (int i=0;i<voxel.nChildren;++i)
   *	    {
   *	      if (visitor.visit(voxel,i))
   *              visitTree(voxel->getChild(i),visitor);
   *
   *          visitor.visited(voxel,i);
   *	    }
   *	  visitor.visited(voxel);
   *    }    
   *  }  
   *
   *  foreach (rootVoxel *voxel)
   *   {
   *     visitor.initialize(voxel);
   *     visitTree_rec(voxel, visitor);
   *  }
   *  \endverbatim
   *
   *  \param nThreads number of threads to use (=glb::num_omp_threads if nThreads<1)
   *  \param visitors A pointer to an array of nThreads visitors (or omp_get_max_threads() 
   *  visitors if nThreads<1): visitors[n] will be assigned to the nth thread.
   *
   *  \tparam V the class of the visitor
   *  \note This version should be used with nThreads>1, when a different visitor needs
   *  to be assigned to each thread (e.g. so that internal variables read/writes do not 
   *  conflict between tasks ...)
   */
  template <class V>
  void visitTree(int nThreads, V visitors[])
  {
    if (nThreads<1) nThreads=glb::num_omp_threads;
    if (nThreads>1)
      {
#pragma omp parallel for num_threads(nThreads) 
	for (int i=0;i<nThreads;++i)
	  {
	    V &visitor = visitors[i];
	    for (unsigned long j=i;j<ROOT_VOXELS_COUNT;j+=nThreads)
	      visitTree_rec(&rootVoxels[j],visitor);
	  }
      }
    else
      {
	for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
	  visitTree_rec(&rootVoxels[i],visitors[0]);
      }
  }


  /** \brief Recursively visit the tree and call visitor.visit(Voxel *voxel) with the 
   *  current voxel as parameter when visiting a voxel for the first time, and 
   *  visitor.visited(Voxel *voxel) when all the descendants of a voxel have been visited.
   *  The descendants of a voxel \a voxel are visited only if visitor.visit(voxel) 
   *  returns true, and the branch is skipped otherwise. Note that if visitor.visit(voxel)
   *  returns false, then visitor.visited(voxel) will not be called. The visitor should 
   *  also implement a 'bool visitor.visit(Voxel *voxel,i)' and 
   *  'void visitor.visited(Voxel *voxel,i)' member that act
   *  similarly, except that they are called before and after visiting the ith child of a 
   *  voxel respectively. The visitor will be used in a way equivalent to:
   *
   *  \verbatim
   *  template <class V>
   *  void visitTree_rec(Voxel *voxel, V &visitor)
   *   {
   *     if (visitor.visit(voxel))
   *     {
   *	  for (int i=0;i<(1L<<NDIM);++i)
   *	    {
   *	      if (visitor.visit(voxel,i))
   *	        visitTree_rec(voxel->getChild(i),visitor);
   *	      visitor.visited(voxel,i);
   *	    }
   *	  visitor.visited(voxel);
   *    }    
   *  }  
   *  \endverbatim
   *
   *  \param visitor the visitor implementing bool V::visit(Voxel*) and 
   *  void V::visited(Voxel*)
   *  \param nThreads number of threads to use
   *  \tparam V the class of the visitor
   *  \warning the private data of the visitor is shared among threads, so if you need 
   *  to update private variable, use visitTree(int nThreads, V visitors[]) instead to
   *  avoid synchronization problems.
   */
  template <class V>
  void visitTree(V & visitor, int nThreads=glb::num_omp_threads)
  {
    if (nThreads>1)
      {
#pragma omp parallel for num_threads(nThreads) 
	for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
	  visitTree_rec(&rootVoxels[i],visitor);
      }
    else
      {
	for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
	  visitTree_rec(&rootVoxels[i],visitor);
      }
  }
  
  /** \brief see visitTree
   */
  template <class V>
  void visitTree(const V & visitor, int nThreads=glb::num_omp_threads)
  {
    if (nThreads>1)
      {
#pragma omp parallel for num_threads(nThreads) 
	for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
	  visitTree_rec(&rootVoxels[i],visitor);
      } 
    else
      {
	for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
	  visitTree_rec(&rootVoxels[i],visitor);
      }
  }
  

  /** \brief return the maximum coordinate of the bounding box along dimension \a d
   */
  double getBBoxMax(int d) const
  {
    return xMax[d];
  }

  /** \brief return the minimum coordinate of the bounding box along dimension \a d
   */
  double getBBoxMin(int d) const
  {
    return xMin[d];
  }

  /** \brief return the sizeof the bounding box along dimension \a d
   */
  double getBBoxSize(int d) const
  {
    return deltaX[d];
  }
  
  /** \brief Compute the integer size of a voxel at level 'level'. 
   *   This size is measured in units of voxels at MAX_LEVEL_FROM_ROOT+1
   *  \param level the level of the voxel
   */
  static ICoord voxelIntLength(int level)
  {
    return ICoord(1)<<(MAX_LEVEL_FROM_ROOT - level + 1);
  }
  
  /** \brief Compute half the integer size of a voxel at level 'level'. 
   *   This size is measured in units of voxels at MAX_LEVEL_FROM_ROOT+1
   *  \param level the level of the voxel
   */
  static ICoord voxelIntHalfLength(int level)
  {
    return (ICoord(1)<<(MAX_LEVEL_FROM_ROOT - level));
  } 
  
  /** \brief Compute the length of a voxel at level 'level' along
   *   a given dimension.
   *  \param level the level of the voxel
   *  \param dim the dimension along which to compute the length
   *  \return The length of the voxel along dimension \a dim
   */
  double getVoxelLength(int level, int dim) const
  {
    return ICoordUnitLen[dim] * voxelIntLength(level);
  }
  
  /** \brief Compute the half length of a voxel at level 'level' 
   *   along a given dimension
   *  \param level the level of the voxel
   *  \param dim the dimension along which to compute the length
   *  \return The half length of the voxel along dimension \a dim
   */  
  double getVoxelHalfLength(int level, int dim) const
  {
    return ICoordUnitLen[dim] * voxelIntHalfLength(level);
  }

  /** \brief returns the integer coordinate along dimension 'dim' of the lower left corner
   *  of a voxel with index 'index' at root level. The coordinate is given within a 
   *  uniform grid at level ROOT_LEVEL.
   */
  static ICoord getRootICoordFromIndex(ICoord index, int dim)
  {
    return ((index>>(INDEX_DEC*dim))&INDEX_MASK)>>MAX_LEVEL_FROM_ROOT_PLUS_ONE;
  }

  /** \brief returns the integer coordinate along dimension 'dim' of the lower left corner
   *  of a voxel with index 'index'. 
   *  The coordinates is given within a uniform grid at level ROOT_LEVEL+levelFromRoot.
   */
  static ICoord getICoordFromIndex(ICoord index, int dim, int levelFromRoot)
  {
    return ((index>>(INDEX_DEC*dim))&INDEX_MASK)>>(MAX_LEVEL_FROM_ROOT_PLUS_ONE - levelFromRoot);
  }

  /** \brief returns the integer coordinates  of the lower left corner
   *  of a voxel with index 'index'. 
   *  The coordinates are given within a uniform grid at level ROOT_LEVEL+levelFromRoot.
   */
  static void getICoordsFromIndex(ICoord out[NDIM], ICoord index, int levelFromRoot=0)
  {
    const ICoord dec=(MAX_LEVEL_FROM_ROOT_PLUS_ONE - levelFromRoot);
    for (int i=0;i<NDIM;++i)
      {
	out[i]=(index&INDEX_MASK)>>dec;
	index>>=INDEX_DEC;
      }
  }

private:

  //! This will init the memory pools and set correct index values for the root voxels.
  //! Should be called by constructor only.
  void initVoxels()
  {
    for (int i=0;i<MAX_TH;++i)
      {
	std::stringstream ss;
	ss<<"voxelGroup@"<<i;
	voxelGroupPool[i].setElementName(ss.str());	
	// As we will allocate voxels by groups of CHILDREN_COUNT, we need to ensure
	// that they will be contiguous in memory (i.e. not split over two pages)
	// voxelPool[i].setGranularity(CHILDREN_COUNT);  
      }
    /*
    std::stringstream ss;
    ss<<"voxelGroup@"<<0;
    voxelGroupPool.setElementName(ss.str());	
    */

    long x[NDIM+1];
    for (int i=0;i<=NDIM;++i) x[i]=0;
    for (unsigned long i=0;i<ROOT_VOXELS_COUNT;++i)
      {
	rootVoxels[i].setEmpty();

	ICoord index = (x[0]*ROOT_VOXEL_ILEN)+ROOT_VOXEL_IHALFLEN;	
	for (int j=1;j<NDIM;j++)
	  index += ((x[j]*ROOT_VOXEL_ILEN)+ROOT_VOXEL_IHALFLEN)<<(INDEX_DEC*j);

	rootVoxels[i].setIndex(index);
	    
	x[0]++;
	if (x[0]>=N_ROOT_PER_DIM) 
	  {
	    x[0]=0;x[1]++;
	    int u=1;
	    while ((x[u]>=N_ROOT_PER_DIM) && (u<NDIM))
	      {x[u]=0;u++;x[u]++;}
	  }
      }    
  }
  
  //! Allocate a new group of (1<<NDIM) voxel at level 'level'
  void popVoxelGroup(Voxel **voxelGroup, int level)
  {
    int th = omp_get_thread_num();
    VoxelGroup *vg;
    voxelGroupPool[th].pop(&vg);
    //voxelGroupPool.pop(&vg);
    (*voxelGroup) = vg->data;
    voxelGroup[0]->poolId=th;
    nLeaves[th] += ((1L<<NDIM)-1);
    //nVoxels += (1L<<NDIM);
    if (th==0)
      {
	verticesAreAssigned=false;
	segmentsAreAssigned=false;
      }
  }

  //! Recycle (delete) a group of (1<<NDIM) voxel at level 'level'
  //! (*voxelGroup) is set to NULL
  //! THIS IS NOT THREAD SAFE !!!!
  void recycleVoxelGroup(Voxel **voxelGroup, int level)
  {
    int th=voxelGroup[0]->poolId;
    VoxelGroup *vg = reinterpret_cast<VoxelGroup*>(*voxelGroup);
    voxelGroupPool[th].recycle(vg);
    //voxelGroupPool.recycle(vg);
    (*voxelGroup)=NULL;
    nLeaves[th] -= ((1L<<NDIM)-1);
    //nVoxels -= (1L<<NDIM);
    if (th==0)
      {
	verticesAreAssigned=false;
	segmentsAreAssigned=false;
      }
  }

private:
  // Internally, we use integer coordinates (iCoords) for the voxels.
  // iCoords are measured in units of half the smallest voxel size,
  // that is the size of a voxel at level (MAX_LEVEL_FROM_ROOT+1).
  // A voxel at level L has integer size iSize=1<<(MAX_LEVEL_FROM_ROOT-L+1) 

  FILE *fl;
 
  void openFile(const char *name)
  {
    static int count =0;
    char fname[255];
    sprintf(fname,"%s_%6.6d.dat",name,count++);
    fl=fopen(fname,"w");
  }

  void closeFile()
  {
    fclose(fl);
  }
  
  /** 
   *  \return the index in the root array of the root voxel with id 'index'
   */
  static ICoord getRootArrayIndex(ICoord index)
  {
    ICoord result=(index&INDEX_MASK)>>MAX_LEVEL_FROM_ROOT_PLUS_ONE;
    ICoord dec = 0;
   
    for (int i=1;i<NDIM;++i)
      {
	index>>=INDEX_DEC;   
	dec+=ROOT_LEVEL;
	result |= ((index&INDEX_MASK)>>MAX_LEVEL_FROM_ROOT_PLUS_ONE)<<dec;	
      }
    
    return result;
  }
   
  /*
  void clearRec(Voxel *voxel)
  {
    if (!voxel->isLeaf())
      {
	for (int i=0;i<(1L<<NDIM);++i)
	  clearRec(voxel->getChild(i));	  

	recycleVoxelGroup(&(voxel->child),voxel->getLevel()+1);

      }
    //voxelGroupPool[voxel->getLevel()].recycle(&voxel); 
  }
  */
  
  /*
  // bIBox[0] is the integer coordinate of vertex 0 of the bounding box
  // bIBox[1] is the integer coordinate of the corner opposite to vertex 0
  // voxIDim[0..NDIM-1] is the integer coordinate of the 0 vertex of the voxel
  // voxIDim[NDIM] is half the integer length of the current voxel
  template <class T, class T2, class OutputIterator>
  static void getLeavesBoxOverlap_rec(const T (&bIBox)[2][NDIM],
				      T2 (&voxIDim)[NDIM+1], 
				      Voxel *voxel, OutputIterator out)
  {    
    if (voxel->isLeaf())
      {
	(*out)=voxel;
	++out;
      }
    else
      {
	for (int i=0;i<(1L<<NDIM);++i)
	  {
	    bool inBound=true;
	    //bool allIn=true;

	    // update voxIDim and check bounding box
	    for (ICoord j=0;j<NDIM;++j)
	      {
		voxIDim[j]+= ((i&(1L<<j))>>j)*voxIDim[NDIM];
		//if (i&(1L<<j)) voxIDim[j]+=voxIDim[NDIM];
						  
		if ((voxIDim[j]>bIBox[1][j])||
		    (voxIDim[j]+voxIDim[NDIM]<bIBox[0][j])) 
		  inBound=false;// out of bound!	
	      }

	    if (inBound) 
	      {
		voxIDim[NDIM]>>=1;
		getLeavesBoxOverlap_rec(bIBox,voxIDim,voxel->getChild(i),out);
		voxIDim[NDIM]<<=1;
	      }

	    // undo what has been done ...
	    for (ICoord j=0;j<NDIM;++j)
	      {	
		voxIDim[j]-= ((i&(1L<<j))>>j)*voxIDim[NDIM];
		//if (i&(1L<<j)) voxIDim[j]-=voxIDim[NDIM];		
	      }	    
	  }
      }
  }
*/

  template <class V>
  void visitTree_rec(Voxel *voxel, V &visitor)
  {
    visitor.initialize(voxel);
    visitTree_recHandler(voxel,visitor);
  }

  template <class V>
  void visitTree_recHandler(Voxel *voxel, V &visitor)
  {
    if (visitor.visit(voxel))
      {
	
	for (int i=0;i<(1L<<NDIM);++i)
	  {
	    if (visitor.visit(voxel,i))
	      visitTree_recHandler(voxel->getChild(i),visitor);
	    visitor.visited(voxel,i);
	  }
	
	visitor.visited(voxel);
      }    
  }  
  
  template <class M, class SH>
  void refineOverSimplex(const M *mesh, const SH s, int level)
  {
    double bBox[2][NDIM];
    mesh->template computeBoundingBox<SH>(s,bBox);
    internal::localAmrGridOverlapVisitor::ThreadedRefineT<MyType> 
      visitor(this,level);		      
    visitBBoxOverlap(bBox,visitor);
  }
  /*
  template <class M, class SH, class V>
  void setRootLevelOverSimplex(const M *mesh, const SH s, int level, 
			   V& rootArr)
  {
    double bBox[2][NDIM];
    mesh->template computeBoundingBox<SH>(s,bBox);
    internal::localAmrGridVisitor::setRootLevelT<MyType> 
      visitor(this,level,rootArr));		      
    visitBBoxOverlap(bBox,visitor);
  }
  */

protected:  
  friend class internal::localAmrGridVisitor::ClearVoxelsT<MyType>;
  friend class internal::localAmrGridVisitor::ToAsciiT<MyType>;
  
  GeometricProperties *geometry;
  //MpiCommunication *mpiCom;

  std::vector<double> xMin;
  std::vector<double> xMax;
  std::vector<double> deltaX;
  std::vector<double> deltaX_inv;

#ifdef HAVE_BOOST
#ifdef HAVE_GMP
  typedef boost::multiprecision::mpf_float mpfloat;
  std::vector<mpfloat> mp_xMin;
  std::vector<mpfloat> mp_xMax;
  std::vector<mpfloat> mp_deltaX;
  std::vector<mpfloat> mp_deltaX_inv;
  mpfloat mp_BBOX_ILEN;
#endif
#endif

  std::vector<ICoord> rootIStride;
  double deltaX_max;
  double deltaX_min;
  std::vector<double> ICoordUnitLen;  
  std::vector<double> epsilon;  
  std::vector<double> halfEpsilon;  
  std::vector<double> halfEpsilonILen;

  unsigned long nLeaves[MAX_TH];
  //unsigned long nVoxels;
  unsigned long nUniqueVertices;
  unsigned long nUniqueSegments;
  bool verticesAreAssigned;
  bool segmentsAreAssigned;

  double voxelVolume[MAX_LEVEL_FROM_ROOT+1];
  double voxelInverseVolume[MAX_LEVEL_FROM_ROOT+1];
  
  //std::vector<Voxel> rootVoxels;
  Voxel rootVoxels[ROOT_VOXELS_COUNT];
  
  struct VoxelGroup
  {
    Voxel data[CHILDREN_COUNT];
  };

  typedef MemoryPoolT<VoxelGroup> VoxelGroupPool;  
  VoxelGroupPool voxelGroupPool[MAX_TH];
};

/** \}*/
#include "../internal/namespace.footer"
#endif
