#ifndef __SOLVER_TOOLS_REFINE_INTERNAL_HXX__
#define __SOLVER_TOOLS_REFINE_INTERNAL_HXX__

#include <limits>

#include "../../internal/namespace.header"

namespace slv { 
  namespace refine { 
    namespace internal {
      
      template <int NDIM, int NDIM_W>
      double computeInvariant(const double dz[2][NDIM_W])
      {
	double result=0;
	for (int i=0;i<NDIM;++i)
	  result+=dz[0][i]*dz[1][NDIM+i] - dz[1][i]*dz[0][NDIM+i];
	return result;
      }

      template <int NDIM, int NDIM_W, typename T>
      double computeOrder2Invariant(const T *v0, const T *v1, const T *v2,
				    const T *t0, const T *t1, const T *t2)
      {
	double result=0;
	T base[2][NDIM_W];
	const T *p[4][3]={{v0,t0,t1},{v1,t2,t0},{v2,t1,t2},{t0,t2,t1}};

	for (int i=0;i<4;++i)
	  {
	    for (int j=0;j<2;++j)
	      for (int k=0;k<NDIM_W;++k)
		base[j][k] = p[i][j+1][k] - p[i][0][k];
	    //geometry->template getBaseVectors<T,T,2,NDIM_W>(p[i],base);
	    result+=computeInvariant<NDIM,NDIM_W>(base);
	  }

	return result;
      }

      // This is just an approximation, recursive version not yet implemented !
      // use with MAX_LEVEL=1 only.
      template <int MAX_LEVEL, class M, class F>
      int findSplitSegmentHelper_REC_OLD(M *mesh, 
					 std::vector<typename M::Simplex>
					 simplices[M::Simplex::NSEG][2],
					 F &functor,
					 double threshold, 
					 int level)
      {
	typedef typename M::Simplex Simplex;
	typedef typename M::Vertex Vertex;

	double globalRet=-1;
	double globalInv=threshold;	
	std::vector<Simplex> newSimplices[Simplex::NSEG][2];

	long nSimplices=0;
	for (int k=0;k<Simplex::NSEG;++k)
	  for (int j=0;j<2;++j)
	    nSimplices+=simplices[k][j].size();

	std::vector<Vertex> newVertices(nSimplices*Simplex::NSEG);
	Vertex *newVertex = &newVertices[0];

	for (int k=0;k<Simplex::NSEG;++k)
	  {
	    double splitInvariant=globalInv;
	    for (int j=0;j<2;++j)
	    {
	      std::vector<Simplex> &curSimplices=simplices[k][j];
	      const long sz=curSimplices.size();
	      
	      // store simplices only if we may need to go to next level
	      if ((globalRet<0)&&(level<MAX_LEVEL))
		newSimplices[k][j].reserve(sz*2*Simplex::NSEG);
	      
	      for (long i=0;i<sz;++i,newVertex+=Simplex::NSEG)
		{
		  //long grp = i>>(level-1);
		  //Vertex newVertex[Simplex::NSEG];
		  std::vector<Simplex*> s[Simplex::NSEG];
		  std::vector<Simplex>  splitS[Simplex::NSEG][2];
		  
		  // THIS MUST BE FALSE -> generate only 2 new simplices.
		  static const bool splitNeighbors=false;
		  mesh->simulateAllSimplexSplits
		    (&curSimplices[i],newVertex,s,splitS,false,true,splitNeighbors);

		  // val stores the best invariant obtained by splitting this simplex
		  double val=threshold;		
		  for (int seg=0;seg<Simplex::NSEG;++seg)
		    {
		      double curInvariant=-1;
		      for (long ss=0;ss<splitS[seg][0].size();++ss)
			{		    
			  // Compute the new invariants		   
			  double nv1 = functor(&splitS[seg][0][ss]);		 
			  if (nv1>curInvariant) curInvariant=nv1;

			  double nv2 = functor(&splitS[seg][1][ss]);
			  if (nv2>curInvariant) curInvariant=nv2;		  
			}	
		      
		      if (val>curInvariant)
			val=curInvariant;		      
		    }	
		  			 
		  // Now we only need to store simplices for next level if we actually
		  // go to the next level ...
		  if ((globalRet<0)&&(level<MAX_LEVEL))
		    {
		      for (int k2=0;k2<Simplex::NSEG;++k2)
			{
			  newSimplices[k][j].push_back(splitS[k2][0][0]);
			  newSimplices[k][j].push_back(splitS[k2][1][0]);
			}
		    }
		  
		  if (val<splitInvariant)
		    splitInvariant=val;
		} // i
	    } // j

	    if (splitInvariant<globalInv)
	      {
		globalInv=splitInvariant;
		globalRet=k;
	      }	    
	  }// k

	if ((globalRet<0)&&(level<MAX_LEVEL))
	  {
	    return findSplitSegmentHelper_REC_OLD<MAX_LEVEL>
	      (mesh,newSimplices,functor,threshold,level+1);
	  }
	
	return globalRet;
      }

      template <int NSEG_LEFT>
      struct FindSplitSegmentHelper_REC
      {
	int maxLevel_;
	double threshold_;	
	double value_;

	char getLevel() {return maxLevel_;}
	double getValue() {return value_;}

	FindSplitSegmentHelper_REC(int maxLevel, double threshold,
				   double value=-1)
	{
	  maxLevel_=maxLevel;
	  threshold_=threshold;	  
	  value_=(value<0)?(std::numeric_limits<double>::max()):value;
	}

	template <class M, class F>
	int compute(M *mesh, 
		    typename M::Simplex *s,
		    typename M::Simplex splitSimplices[M::Simplex::NSEG][2],
		    F &functor)
	{
	  typedef typename M::Simplex Simplex;
	  int ret=-1;
	  double val=threshold_;
	  int level=Simplex::NSEG;
	  char segPerm[Simplex::NSEG-1];

	  for (int seg=0;seg<Simplex::NSEG;++seg)
	    {
	      for (int i=0,j=0;j<Simplex::NSEG;++j)
		if (seg!=j) segPerm[i++]=j;
	      
	      get(mesh,s,&splitSimplices[seg][0],functor,segPerm);
	      double newVal=getValue();
	      double newLvl=getLevel();
	      
	      if (newVal < threshold_)
		{
		  if ((newLvl<level)||(val>newVal))
		    {
		      val=newVal;
		      level=newLvl;
		      ret=seg;
		    }
		}
	    }

	  return ret;
	}

	template <class M, class F, int LEVEL=1>
	void get(M *mesh, 
		typename M::Simplex *s,
		typename M::Simplex simplices[2+LEVEL-1],
		F &functor,		       
		char curSegPerm[NSEG_LEFT])
	{
	  static const int N_SIMPLICES = (2+LEVEL-1);

	  typedef typename M::Simplex Simplex;
	  typedef typename M::Vertex Vertex;
	  typedef typename Simplex::SegmentHandle SegmentHandle;
	  
	  // Try splitting every remaining segment
	  for (int seg=0;seg<NSEG_LEFT;++seg)
	    {	    	    
	      // That's the segment we have to split first
	      SegmentHandle handle=s->getSegmentHandle(curSegPerm[seg]);	    	    
	      
	      // find which simplex it belongs to ...
	      for (int sid=0;sid<N_SIMPLICES;sid++)
		{
		  int segId=simplices[sid].findSegmentIndex(handle);	  
		  if (segId>=0) 
		    {
		      // That's the one
		      Simplex splitSimplices[2];
		      Vertex newVertex;
		      //printf("REC (v=%e, ml=%d, cur= %d)\n",value_,maxLevel_,LEVEL);
		      // simulate splitting along this segment
		      mesh->simulateIsolatedSimplexSplit
			(&simplices[sid],segId,newVertex,splitSimplices,false,true);

		      // compute the new invariant
		      double val=0;
		      for (int j=0;j<N_SIMPLICES;++j)
			if (j!=sid)
			  {
			    double nv=functor(&simplices[j]);
			    if (nv>val) val=nv;
			  }
		      for (int i=0;i<2;++i)
			{
			  double nv=functor(&splitSimplices[i]);
			  if (nv>val) val=nv;
			}
		      
		      // check if it is an improvement
		      if (val<threshold_)
			{
			  if ((LEVEL<maxLevel_)||(val<value_))
			    {
			      value_=val;			      
			      maxLevel_=LEVEL;
			    }
			}		    
		      
		      // and continue to next segment if needed
		      if ((value_ >= threshold_)&&(LEVEL<maxLevel_))
			{
			  // setup new set of simplices by adding the split ones
			  Simplex newSimplices[N_SIMPLICES+1];
			  char newSegPerm[NSEG_LEFT-1];
			  // the new set of remaining segments to try
			  for (int i=0,j=0;j<NSEG_LEFT;++j)
			    if (j!=seg) newSegPerm[i++]=curSegPerm[j];

			  std::copy_n(simplices,N_SIMPLICES,newSimplices);
			  newSimplices[sid]=splitSimplices[0];
			  newSimplices[N_SIMPLICES]=splitSimplices[1];
		      
			  // recurse !
			  FindSplitSegmentHelper_REC<NSEG_LEFT-1> 
			    next(maxLevel_,threshold_,value_);
		      			  
			  next.template get<M,F,LEVEL+1>
			    (mesh,s,newSimplices,functor,newSegPerm);

			  // and check if recursing improved the situation
			  if ((next.value_ < value_)||
			      (next.maxLevel_<maxLevel_))
			    {			      
			      value_=next.value_;
			      maxLevel_=next.maxLevel_;
			    }
			}

		      break;
		    }
		}
	    }	
	}
      };

      template <>
      struct FindSplitSegmentHelper_REC<0>
      {
	int maxLevel_;
	double threshold_;
	double value_;

	FindSplitSegmentHelper_REC(int maxLevel, double threshold,
				   double value=-1)
	{
	  maxLevel_=maxLevel;
	  threshold_=threshold;
	  value_=(value<0)?(std::numeric_limits<double>::max()):value;
	}

	template <class M, class F, int LEVEL>
	void get(M *mesh, 
		 typename M::Simplex *s,
		 typename M::Simplex simplices[2+LEVEL],
		 F &functor, 
		 char *curSegPerm)
	{}
      };

      template <class M, class F>
      int findSplitSegmentHelper(M *mesh, typename M::Simplex *s, 
				 F &functor,double threshold, 
				 bool recursive=false)
      {
	typedef typename M::Simplex Simplex;
	typedef typename M::Vertex Vertex;

	Vertex newVertex[Simplex::NSEG];	
	Simplex splitSimplices[Simplex::NSEG][2];
	
	mesh->simulateAllIsolatedSimplexSplits
	  (s,newVertex,splitSimplices,false,true);

	double val=threshold;
	int ret=-1;	   

	// Try finding the segment that improves the invariant best when split
	for (int seg=0;seg<Simplex::NSEG;++seg)
	  {
	    double curInvariant=0;
	    
	    // Compute the new invariants		   
	    double nv1 = functor(&splitSimplices[seg][0]);		 
	    if (nv1>curInvariant) curInvariant=nv1;

	    double nv2 = functor(&splitSimplices[seg][1]);
	    if (nv2>curInvariant) curInvariant=nv2;		  	  
		
	    // And keep the segment that minimizes the max of the new invariants
	    if (val>curInvariant)
	      {
		val=curInvariant;
		ret=seg;		
	      }
	  }

	/*
	if (ret < 0)
	  {
	    std::vector<Simplex> simplices[Simplex::NSEG][2];
	    for (int i=0;i<Simplex::NSEG;++i)
	      {
		simplices[i][0].push_back(splitSimplices[i][0]);
		simplices[i][1].push_back(splitSimplices[i][1]);
	      }
	    ret = 
	      findSplitSegmentHelper_REC_OLD<1>
	      (mesh,simplices,functor,threshold,1);	    	    
	  }
	*/	
	
	if ((ret<0)&&(recursive))
	  {
	    // limit first search to 2 levels of recursion
	    int maxLevel=std::min(Simplex::NSEG-1,2);
	    FindSplitSegmentHelper_REC<Simplex::NSEG-1>
	      findSegmentHelper(maxLevel,threshold);

	    ret=findSegmentHelper.compute(mesh,s,splitSimplices,functor);

	    if ((ret<0)&&(maxLevel<Simplex::NSEG-1))
	      {
		// full search only if needed
		maxLevel=Simplex::NSEG-1;
		FindSplitSegmentHelper_REC<Simplex::NSEG-1> 
		  findSegmentHelper2(maxLevel,threshold);

		ret=findSegmentHelper2.compute(mesh,s,splitSimplices,functor);
	      }
	  }
	

	return ret;
      }

      template <class M, class F>
      int findSplitSegmentHelper_OLD(M *mesh, 
				 typename M::Simplex *s, 
				 F &functor,
				 double threshold)
      {
	typedef typename M::Simplex Simplex;
	typedef typename M::Vertex Vertex;

	Vertex newVertex[Simplex::NSEG];
	std::vector<Simplex*> simplices[Simplex::NSEG];
	std::vector<Simplex>  splitSimplices[Simplex::NSEG][2];
	
	// Last param sets whether we check the simplex new invariant (false) or
	// all the affected neighborhood (true)
	static const bool splitNeighbors=false;
	mesh->simulateAllSimplexSplits
	  (s,newVertex,simplices,splitSimplices,false,true,splitNeighbors);

	double val=threshold;
	int ret=-1;
	   
	// double neiInv[Simplex::NNEI];
	// if (splitNeighbors)
	//   {
	// 	for (int nei=0;nei<Simplex::NNEI;++nei)
	// 	  {
	// 	    Simplex *neighbor=s->getNeighbor(nei);
	// 	    if (neighbor!=NULL)
	// 	      neiInv[nei]=slv::refine::
	// 		poincareInvariantFromTracers<Mesh>(s->getNeighbor(nei),geometry);
	// 	    else neiInv[nei]=-1;
	// 	  }
	//   }

	// Try finding the segment that improves the invariant best when split
	for (int seg=0;seg<Simplex::NSEG;++seg)
	  {
	    double curInvariant=0;
	    for (long i=0;i<splitSimplices[seg][0].size();++i)
	      {		    
		// Compute the new invariants		   
		double nv1 = functor(&splitSimplices[seg][0][i]);		 
		if (nv1>curInvariant) curInvariant=nv1;

		double nv2 = functor(&splitSimplices[seg][1][i]);
		if (nv2>curInvariant) curInvariant=nv2;		  
	      }			
		
	    // Check also the value of the unaffected neighbors
	    // if (splitNeighbors)
	    //   {
	    //     SegmentHandle handle=s->getSegmentHandle(seg);
	    //     for (int nei=0;nei<Simplex::NNEI;++nei)
	    //       {
	    // 	Simplex *neighbor=s->getNeighbor(nei);
	    // 	if (neiInv[nei]>curInvariant) 
	    // 	  {
	    // 	    if ((neighbor->getVertexIndex(handle->getVertex(0))<0)||
	    // 		(neighbor->getVertexIndex(handle->getVertex(1))<0))
	    // 	      curInvariant=neiInv[nei];	
	    // 	  } 
	    //       }	    
	    //   }
		
	    // And keep the segment that minimizes the max of the new invariants
	    if (val>curInvariant)
	      {
		val=curInvariant;
		ret=seg;
		//ret.second=seg;
	      }
	  }
	

	if (ret < 0)
	  {
	    ret = 
	      findSplitSegmentHelper_REC_OLD<1>
	      (mesh,splitSimplices,functor,threshold,1);
	    
	    // ret = 
	    //   findSplitSegmentHelper_REC<Simplex::NSEG-1>
	    //   (mesh,splitSimplices,functor,threshold,1);
	    
	  }
	
	return ret;
      }
  
    }
  }
}

#include "../../internal/namespace.footer"
#endif
