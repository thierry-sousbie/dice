#ifndef __LOCAL_AMR_GRID_DEFAULT_VISITORS_HXX__
#define __LOCAL_AMR_GRID_DEFAULT_VISITORS_HXX__

#include <stdio.h>
#include <stack>

#include "../../internal/namespace.header"

namespace internal {
  namespace localAmrGridVisitor {
  
    template <class AMR>
    class ClearVoxelsT
    {
      typedef typename AMR::Voxel Voxel;    
    public:
      ClearVoxelsT(AMR *amr_):amr(amr_)
      {}

      static void initialize(Voxel *root){}

      static bool visit(Voxel *voxel)
      {      
	return (!voxel->isLeaf());     
      }

      static bool visit(Voxel *voxel,int i) 
      {
	return true;
      }

      void visited(Voxel *voxel) const
      {
	amr->recycleVoxelGroup(&(voxel->child),voxel->getLevel()+1);     
      }

      static void visited(Voxel *voxel, int i) 
      {
      }
    
    private:
      mutable AMR *amr;    
    };

    template <class AMR>
    class AssignVerticesT
    {
      typedef typename AMR::Voxel Voxel;
      static const long NDIM=AMR::NDIM;
    public:
      AssignVerticesT():
	amr(NULL),
	assignSegments(false) 
      {nVert=0;nSeg=0;}
      AssignVerticesT(AMR *amr_, bool assignSegments_):
	amr(amr_),
	assignSegments(assignSegments_)      
      {nVert=0;nSeg=0;}

      static void initialize(Voxel *root){}

      bool visit(Voxel *voxel)
      { 
	if (voxel->isLeaf())
	  {
	    Voxel *nei[1L<<NDIM];
	    int level=voxel->getLevel();

	    voxel->resetFlags();
	 	  
	    for (int vertexId=0;vertexId<(1L<<NDIM);++vertexId)
	      {
		bool vertexIsMine=true;	      

		amr->getVertexNeighbors(voxel,vertexId,nei);

		for (int voxelId=1;voxelId<(1L<<NDIM);++voxelId)
		  {
		    // if ((nei[voxelId]->getLevel() > level) ||
		    // 	((nei[voxelId]->getLevel()==level)&&(nei[voxelId]>voxel)))

		    // A voxel owns a vertex if, among all neighbor voxels sharing it, its 
		    // level is the highest and no neighbor at same level have a 
		    // higher index.
		    if ((nei[voxelId]->getLevel() > level) ||
			((nei[voxelId]->getLevel()==level)&&
			 (nei[voxelId]->getIndex()>voxel->getIndex())))
		      vertexIsMine=false;
		  }

		if (vertexIsMine) 
		  {
		    // if (check) printf("Adding flag for vertex %d\n",i);
		    voxel->setVertexFlag(vertexId);
		    ++nVert;
		  }

		if (assignSegments)
		  {
		    for (int d=0;d<NDIM;++d)
		      {		  
			if ((vertexId&(1<<d))==0)
			  {
			    bool segmentIsMine=true;
			    for (int voxelId=1;voxelId<(1L<<NDIM);++voxelId)
			      {
				// (voxelId&(1<<d)) is true if the voxel was shifted
				// along dimension 'd' -> we need those voxels that
				// were not shifted along 'd'.
				// if (
				//     ((voxelId&(1<<d))==0) &&
				//     ((nei[voxelId]->getLevel() > level) ||
				//      ((nei[voxelId]->getLevel()==level)&&
				//       (nei[voxelId]>voxel)) )
				//     )
				if (
				    ((voxelId&(1<<d))==0) &&
				    ((nei[voxelId]->getLevel() > level) ||
				     ((nei[voxelId]->getLevel()==level)&&
				      (nei[voxelId]->getIndex()>voxel->getIndex())) )
				    )
				  segmentIsMine=false;			    
			      }
			    if (segmentIsMine) 
			      {			 
				voxel->setSegmentFlag(vertexId,d);
				++nSeg;
			      }
			  }
		      }
		  }

		/*
		  if (mine)
		  {
		  double corner[NDIM];
		  double center[NDIM];
		  amr->index2CenterCoords(voxel->getIndex(),center);
		  amr->index2CornerCoords(voxel->getIndex(),voxel->getLevel(),corner,i);
		  fprintf(fl,"%lg %lg %lg %lg\n",corner[0],corner[1],center[0],center[1]);
		  }
		*/
	    
	      }
	    return false;
	  }
	return true;      
      }
    
      static bool visit(Voxel *voxel,int i) 
      {
	return true;
      }
    
      static void visited(Voxel *voxel) 
      {
      }
    
      static void visited(Voxel *voxel, int i) 
      {
      }

      void init(AMR *amr_, bool assignSegments_) 
      {
	amr=amr_;assignSegments=assignSegments_;
	nVert=0;nSeg=0;
      }
      unsigned long getNVertices() const {return nVert;}
      unsigned long getNSegments() const {return nSeg;}
    private:
      mutable AMR *amr;
      mutable bool assignSegments;
      unsigned long nVert;
      unsigned long nSeg;
    };

    namespace internal {
      template <int NDIM> struct GetBoundaryFlagsHelper;      
    }

    // Here pass the level of neighbors in order -> 
    // construct NDIMS boundary flags from that !
    // Boundary flags work on 32 bits, with the left/right  boundary of dimension D at 
    // 1<<D<<16 and 1<<D respectively.
    template <class AMR>
    class AssignVertices_FastT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM=AMR::NDIM;
      static const unsigned int MASK = (1<<NDIM)-1;
    public:
      AssignVertices_FastT():
	amr(NULL),
	assignSegments(false) 
      {nVert=0;nSeg=0;}
      AssignVertices_FastT(AMR *amr_,bool assignSegments_,
			   const Voxel *rootVoxels,
			   const unsigned char *rootVoxelsLevel):
	amr(amr_),	
	assignSegments(assignSegments_),
	slowVisitor(amr_,assignSegments_),
	ref(rootVoxels),
	rootLevel(rootVoxelsLevel)
      {nVert=0;nSeg=0;}

      void initialize(Voxel *root)
      {
	boundaryFlagsAny=0;
	locationFlags= MASK | (MASK<<16);//(~0);
	globalBoundaryFlags = getBoundaryFlags(root,boundaryFlags);
	boundaryFlagsAny|=boundaryFlags[0];
	for (int i=1;i<NDIM;++i) 
	  boundaryFlagsAny|=boundaryFlags[i];
      }

      bool visit(Voxel *voxel)
      { 
	if (voxel->isLeaf())
	  {
	    voxel->resetFlags();	    	    
	    
	    // A bit is set in locationFlags only if the voxel lies on the corresponding 
	    // boundary, while boundaryFlags[d] bits are set only if root neighbors level
	    // conflicts with the basic vertex/segment attribution procedure.
	    unsigned int genericConfiguration=0;
	    if (locationFlags)
	      {
		// Configuration is always considered generic on the boundaries
		// of the domain. We have to do that because our basic configuration
		// does not correspond to the generic one on boundaries due to the neighbor
		// indexing ... (this should be OK in non periodic case though)
		if (globalBoundaryFlags&locationFlags) genericConfiguration=true;
		else if (boundaryFlagsAny&locationFlags)
		  {
		    // We are in a generic configuration if the location flags intersect 
		    // the boundary flags
		    int nBits=0;
		    for (unsigned int tmp=locationFlags;(tmp!=0);tmp &= (tmp-1))
		      genericConfiguration |= (boundaryFlags[nBits++] & locationFlags);
		    
		    
		    // If we want to treat only the truely generic case, then we still 
		    // need to filter.
		    if (genericConfiguration)
		      {			
			genericConfiguration=(boundaryFlags[0] & locationFlags);
			for (int i=1;i<nBits;++i)
			  {
			    int n=i+1;
			    for (unsigned int tmp=(boundaryFlags[i] & locationFlags);
				 (tmp!=0);
				 tmp &= (tmp-1))
			      --n;
			    // This means we are on at least nBits boundaries, and 
			    // they are all flagged in the boundaryFlags of order nBits 
			    // => we are indeed in a generic configuration
			    // Note that 'n' may be negative if the voxel is a root
			    if (n<=0) genericConfiguration=true;			
			  }
		      }	
		    	
		  }		
	      }
	    
	    // Basic configuration (neighbors have non conflicting level !)
	    if (!genericConfiguration) 
	      {
		// => a basic voxel always owns it's 0th vertex
		voxel->setVertexFlag(0);
		++nVert;

		if (assignSegments) 
		  {
		    // => and the 3 segments departing from it
		    for (int d=0;d<3;++d) 
		      voxel->setSegmentFlag(0,d);
		    nSeg+=3;
		  }
	      }
	    else 
	      {
		//generic configuration (neighbors have conflicting levels)
		unsigned long nv=slowVisitor.getNVertices();
		unsigned long ns=slowVisitor.getNSegments();

		slowVisitor.visit(voxel);

		nVert +=slowVisitor.getNVertices() - nv;
		nSeg  +=slowVisitor.getNSegments() - ns;
	      }

	    return false;
	  }	
	
	locationStack.push(locationFlags);
	return true;      
      }
    
      bool visit(Voxel *voxel,int i) 
      {	
	// remove a location boundary flag if we went in the opposite direction at 
	// least once
	locationFlags = locationStack.top() & ((i | (((~i)&MASK)<<16)));
	return true;
      }
    
      void visited(Voxel *voxel) 
      {
	locationStack.pop();
      }
    
      static void visited(Voxel *voxel, int i) {}
      
      void init(AMR *amr_, bool assignSegments_,
		const Voxel *rootVoxels,
		const unsigned char *rootVoxelsLevel)
      {
	amr=amr_;
	assignSegments=assignSegments_;
	nVert=0;nSeg=0;
	slowVisitor.init(amr_,assignSegments_);
	ref=rootVoxels;
	rootLevel=rootVoxelsLevel;
      }
      
      unsigned long getNVertices() const {return nVert;}
      unsigned long getNSegments() const {return nSeg;}

    protected:
 
      unsigned int getBoundaryFlags(const Voxel *rootVoxel, unsigned int flags[NDIM])
      {
	static const long rootGridSize = (1<<AMR::ROOT_LEVEL);	
	static const long rootGridSizeMinusOne = rootGridSize-1;	

	ICoord iCoord[NDIM];
	long delta[NDIM][2];
	unsigned int globalFlags=0;
	
	//long iStride[NDIM+1];
	amr->getICoordsFromIndex(iCoord,rootVoxel->getIndex());
		
	long stride=1;
	for (int i=0;i<NDIM;++i)
	  {
	    delta[i][0]=-stride;//-iStride[i];
	    delta[i][1]=stride;//iStride[i];
	    stride*=rootGridSize;
	  }

	if (AMR::PERIODIC_BOUNDARIES)
	  {
	    for (int i=0;i<NDIM;++i)
	      {
		if (iCoord[i]==0) 
		  {delta[i][0]*=(-rootGridSizeMinusOne);globalFlags|=((1<<i)<<16);}
		else if (iCoord[i]==(rootGridSizeMinusOne))
		  {delta[i][1]*=(-rootGridSizeMinusOne);globalFlags|=(1<<i);}
	      }
	  }
	else
	  {
	    for (int i=0;i<NDIM;++i)
	      {
		if (iCoord[i]==0) 
		  {delta[i][0]=0;globalFlags|=((1<<i)<<16);}
		else if (iCoord[i]==(rootGridSizeMinusOne))
		  {delta[i][1]=0;globalFlags|=(1<<i);}
	      }
	  }
	
	long index = std::distance(ref,rootVoxel);
	internal::GetBoundaryFlagsHelper<NDIM>::
	  compute(rootLevel+index,delta,flags);

	return globalFlags;
      }

    private:
      mutable AMR *amr;
      mutable bool assignSegments;
      unsigned long nVert;
      unsigned long nSeg;

      std::stack<unsigned int> locationStack;
      unsigned int locationFlags;
      unsigned int boundaryFlags[NDIM];
      unsigned int boundaryFlagsAny;
      unsigned int globalBoundaryFlags;
     
      AssignVerticesT<AMR> slowVisitor;
      const Voxel *ref;
      const unsigned char *rootLevel;
    };

    template <class AMR>
    class ToAsciiT
    {
      typedef typename AMR::Voxel Voxel;
      static const long NDIM=AMR::NDIM;
    public:
      ToAsciiT(AMR *amr_,FILE *f):amr(amr_),f(f)
      {}

      static void initialize(Voxel *root){}

      bool visit(Voxel *voxel) const
      { 
	if (voxel->isLeaf())
	  {
	    double coords[NDIM];
	    double dx[NDIM];
	    for (int i=0;i<NDIM;++i)
	      dx[i]=amr->voxelIntHalfLength(voxel->getLevel())*amr->ICoordUnitLen[i];

	    amr->index2CenterCoords(voxel->getIndex(),coords);    

	    fprintf(f,"%lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		    coords[1]-dx[1],coords[1]-dx[1],
		    coords[1]+dx[1],coords[1]+dx[1],
		    coords[0]-dx[0],coords[0]+dx[0],
		    coords[0]+dx[0],coords[0]-dx[0],
		    (double)voxel->data);
	    return false;
	  }
	return true;      
      }

      static bool visit(Voxel *voxel,int i) 
      {
	return true;
      }

      static void visited(Voxel *voxel) 
      {
      }

      static void visited(Voxel *voxel, int i) 
      {
      }
    
    private:
      mutable AMR *amr;
      mutable FILE *f;
    };
   
    namespace internal {
      namespace {
	template <int NDIM,int D>
	struct GetBoundaryFlagsHelperRec
	{	
	  static void compute(unsigned char refLevel,
			      const unsigned char *curLevel,
			      const long delta[NDIM][2],
			      unsigned int *flags,
			      unsigned int tag)
	  {
	    GetBoundaryFlagsHelperRec<NDIM,D+1>::compute
	      (refLevel,curLevel+delta[D][0],delta,flags+1,tag|((1<<D)<<16));
	    
	    GetBoundaryFlagsHelperRec<NDIM,D+1>::compute
	      (refLevel,curLevel,delta,flags,tag);	    
	    
	    GetBoundaryFlagsHelperRec<NDIM,D+1>::compute
	      (refLevel,curLevel+delta[D][1],delta,flags+1,tag|(1<<(D)));
	  }
	};

	template <int NDIM>
	struct GetBoundaryFlagsHelperRec<NDIM,NDIM>
	{	
	  static void compute(unsigned char refLevel,
			      const unsigned char *curLevel,
			      const long delta[NDIM][2],
			      unsigned int *flags,
			      unsigned int tag)
	  {
	    if ( (*curLevel) != refLevel )
	      (*flags) |= tag;		
	  }
	};
      } // unnamed namespace 

      template <int NDIM>
      struct GetBoundaryFlagsHelper
      {	
	static const int D = 0;
	static void compute(const unsigned char *refLevel,
			    const long delta[NDIM][2], 
			    unsigned int *flags)
	{
	  for (int i=0;i<NDIM;++i) flags[i]=0;
	  
	  GetBoundaryFlagsHelperRec<NDIM,0>::compute
	    (*refLevel,refLevel,delta,flags-1,0);	  
	}
      };    
      
    } // internal(2) namespace
  } // localAmrGridVisitor
} // internal namesapce

#include "../../internal/namespace.footer"

#endif
