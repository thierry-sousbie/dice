#ifndef __SOLVER_TOOLS_REFINE_HXX__
#define __SOLVER_TOOLS_REFINE_HXX__

#include "../dice_globals.hxx"
#include "../mesh/cellData/cellDataFunctors_interface.hxx"
#include "../solver/internal/solverTools_refine_internal.hxx"

/**
 * @file 
 * @brief A set of helper functions that may be usefull for refining the mesh in the 
 * context of solving vlasov equations
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup SOLVER
 *   \{
 */

namespace slv {
  namespace refine {    
    
    template <class M, class F>
    int findSplitSegment(M *mesh,
			 typename M::Simplex *simplex, 
			 F &functor, double threshold,
			 bool recursive=false)
    {
      return internal::findSplitSegmentHelper(mesh,simplex,functor,threshold,recursive);
    }
   
    /** 
     * \brief Identifies the longest segment and return its index and square length
     */
    template <class M>
    typename M::CheckRefineReturnType 
    length2(typename M::Simplex *s, const typename M::GeometricProperties *geometry, 
	    double threshold=0, bool computeSegment=true)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::Simplex Simplex;
      typedef typename M::Coord Coord;
    
      int segIndex=0;
      SegmentHandle ref = s->getSegmentHandle(segIndex);
      double score=geometry->template 
	distance2<Coord,M::NDIM_W>(ref->getVertex(0)->getCoordsConstPtr(),
				   ref->getVertex(1)->getCoordsConstPtr());
    
      for (int i=1;i<Simplex::NSEG;i++)
	{
	  SegmentHandle seg=s->getSegmentHandle(i);
	  double v=geometry->template 
	    distance2<Coord,M::NDIM_W>(seg->getVertex(0)->getCoordsConstPtr(),
				       seg->getVertex(1)->getCoordsConstPtr());
	  if (v==score)
	    {
	      ref = s->getSegmentHandle(segIndex);
	      if ( (*ref) < (*seg) )
		segIndex=i;	
	    }
	  else if (v>score) {score=v;segIndex=i;}
	}
      
      return CheckRefineReturnType(score,segIndex);
    }

    /** 
     * \brief Identifies the longest segment in lagrangian coordinates 
     * (i.e. using initCoords) and return its index and square length
     */
    template <class M>
    typename M::CheckRefineReturnType 
    lagrangianLength2(typename M::Simplex *s, 
		      const typename M::GeometricProperties *geometry, 
		      double threshold=0, bool computeSegment=true)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::Simplex Simplex;
      typedef typename M::Coord Coord;
    
      int segIndex=0;
      SegmentHandle ref = s->getSegmentHandle(segIndex);
      double score=geometry->template 
	distance2<Coord,M::NDIM_W>(ref->getVertex(0)->initCoords.getPointer(),
				   ref->getVertex(1)->initCoords.getPointer());
    
      for (int i=1;i<Simplex::NSEG;i++)
	{
	  SegmentHandle seg=s->getSegmentHandle(i);
	  double v=geometry->template 
	    distance2<Coord,M::NDIM_W>(seg->getVertex(0)->initCoords.getPointer(),
				       seg->getVertex(1)->initCoords.getPointer());
	  if (v==score)
	    {
	      ref = s->getSegmentHandle(segIndex);
	      if ( (*ref) < (*seg) )
		segIndex=i;	
	    }
	  else if (v>score) {score=v;segIndex=i;}
	}
      
      return CheckRefineReturnType(score,segIndex);
    }

    
    template <class M>
    std::pair<double,int>
    poincareInvariant(typename M::Simplex* simplex,
		      const typename M::GeometricProperties *geometry)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
      static const int NVERT = Simplex::NVERT;

      double result=0;
      int index=0;

      const Coord *p[NVERT];
      simplex->getVerticesCoordsConstPtr(p);
      Coord base[NDIM][NDIM_W];
      geometry->template getBaseVectors<Coord,Coord,NDIM,NDIM_W>(p,base);
      
      if (NDIM==2)
	{
	  result=fabs(internal::computeInvariant<NDIM,NDIM_W>(base));	
	}
      else if (NDIM==3)
	{
	  Coord dz[2][NDIM_W];

	  std::copy_n(&base[0][0],NDIM_W,&dz[0][0]);
	  std::copy_n(&base[1][0],NDIM_W,&dz[1][0]);
	  result=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));	
	  index=3;

	  //std::copy_n(&base[0][0],NDIM_W,&dz[0][0]);
	  std::copy_n(&base[2][0],NDIM_W,&dz[1][0]);
	  double pi=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));	
	  if (pi>result) {result=pi;index=2;}

	  std::copy_n(&base[1][0],NDIM_W,&dz[0][0]);
	  //std::copy_n(&base[2][0],NDIM_W,&dz[1][0]);
	  pi=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	  if (pi>result) {result=pi;index=1;}

	  for (int i=0;i<NDIM_W;++i)
	    {
	      dz[0][i]=base[2][i]-base[0][i];
	      dz[1][i]=base[1][i]-base[0][i];
	    }
	  pi=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	  if (pi>result) {result=pi;index=0;}
	}
      else
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("poincare invariant no implemented for NDIM>3.\n");
	  exit(-1);
	}

      return std::make_pair(result,index);     
    }

    template <class M>
    double
    poincareInvariantWithSegTracers(typename M::Simplex* simplex,
				    const typename M::GeometricProperties *geometry)
    {      
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      typedef typename Simplex::SegTracers::Type Tracer;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
      static const int NVERT = Simplex::NVERT;
      static const int NSEG = Simplex::NSEG;

      //const Coord *v[NVERT];
      Coord v[NVERT][NDIM_W];
      Coord t[NSEG*NDIM_W];
      double result=0;

      const Coord *ref = simplex->getVertex(0)->getCoordsConstPtr();
      for (int i=0;i<NSEG;++i)
	geometry->getVector(ref,simplex->segTracers.getPointer()+i*NDIM_W,t+i*NDIM_W);
      for (int i=0;i<NVERT;++i)
	geometry->getVector(ref,simplex->getVertex(i)->getCoordsConstPtr(),v[i]);

      //std::copy_n(simplex->tracer.getPointer(),NSEG*NDIM_W,t);
      //simplex->getVerticesCoordsConstPtr(v);
  
      if (NDIM==2)
	{
	  result = fabs(internal::computeOrder2Invariant<NDIM,NDIM_W>
			(v[0],v[1],v[2],
			 t + NDIM_W*Simplex::template findSegmentIndex<0,1>(),
			 t + NDIM_W*Simplex::template findSegmentIndex<0,2>(),
			 t + NDIM_W*Simplex::template findSegmentIndex<1,2>()));	  
	}
      else if (NDIM==3)
	{	
	  result = fabs(internal::computeOrder2Invariant<NDIM,NDIM_W>
			(v[0],v[1],v[2],
			 t + NDIM_W*Simplex::template findSegmentIndex<0,1>(),
			 t + NDIM_W*Simplex::template findSegmentIndex<0,2>(),
			 t + NDIM_W*Simplex::template findSegmentIndex<1,2>()));
	  
	  double tmp;
	  tmp = fabs(internal::computeOrder2Invariant<NDIM,NDIM_W>
		     (v[0],v[3],v[1],
		      t + NDIM_W*Simplex::template findSegmentIndex<0,3>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<0,1>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<3,1>()));
	  if (tmp>result) result=tmp;

	  tmp = fabs(internal::computeOrder2Invariant<NDIM,NDIM_W>
		     (v[0],v[2],v[3],
		      t + NDIM_W*Simplex::template findSegmentIndex<0,2>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<0,3>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<2,3>()));
	  if (tmp>result) result=tmp;

	  tmp = fabs(internal::computeOrder2Invariant<NDIM,NDIM_W>
		     (v[1],v[3],v[2],
		      t + NDIM_W*Simplex::template findSegmentIndex<1,3>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<1,2>(),
		      t + NDIM_W*Simplex::template findSegmentIndex<3,2>()));
	  if (tmp>result) result=tmp;
	}
      else
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("poincare invariant no implemented for NDIM>3.\n");
	  exit(-1);
	}

      return result;
    }

    template <class M>
    std::pair<double,int>
    poincareInvariantWithSegTracers_order1(typename M::Simplex* simplex,
					   //const typename M::SegmentHandle &seg, int sid,
					   const typename M::GeometricProperties *geometry, 
					   bool average=false)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
      
      Coord dz[2][NDIM_W];
      Coord *coords[3];
      double result=0;
      double avg=0;
      int index=0;

      for (int sid=0;sid<Simplex::NSEG;++sid)
	{
	  //int sid = seg->getSimplex()->findSegmentIndex(seg);
	  SegmentHandle seg=simplex->getSegmentHandle(sid);
	  coords[0]=seg->getVertex(0)->getCoordsPtr();
	  coords[1]=seg->getVertex(1)->getCoordsPtr();
	  
	  for (int sid2=0;sid2<Simplex::NSEG;++sid2)
	    {	      
	      if (sid2!=sid)
		{
		  coords[2]=&simplex->segTracers.getPointer()[sid2*NDIM_W];
		  // poincaré invariant computed from vertex/opp facet barycenter/tracer
		  geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
		  double pi=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
		  if (pi>result) 
		    {
		      result=pi;
		      index=sid;
		    }
		  avg+=pi;
		}
	    }
	}
      if (average) return std::make_pair(avg/(2*Simplex::NSEG),index);
      else return std::make_pair(result,index);
    }

    template <class M>
    typename M::CheckRefineReturnType 
    poincareInvariantFromNeighbors(const typename M::Simplex::Neighborhood &nb, 
		      const typename M::GeometricProperties *geometry, 
		      double threshold=0, bool computeSegment=true)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
  
      SegmentHandle segment;
      double score=0;
      int kMax=0;    
  
      Coord midPoint[NDIM_W];
      Coord dz[2][NDIM_W];
      Coord *coords[3];
    
      for (int k=0;k<Simplex::NNEI;++k)
	{
	  if ((M::BOUNDARY_TYPE !=  BoundaryType::PERIODIC)&&(nb.neiV[k]==NULL)) continue;
	
	  SegmentHandle sh[Simplex::Facet::NSEG];	  

	  nb.fh[k]->getSegmentHandles(sh);
	  coords[0]=nb.curV[k]->getCoordsPtr();
	  coords[1]=nb.neiV[k]->getCoordsPtr();	
	  coords[2]=midPoint;
	
	  for (int j=0;j<Simplex::Facet::NSEG;++j)
	    {
	      geometry->midPointCoords(sh[j]->getVertex(0)->getCoordsPtr(),
				       sh[j]->getVertex(1)->getCoordsPtr(),
				       midPoint);
   
	      geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
	    
	      double invariant = fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	      /*
	      double invariant = 0;
	      for (int i=0;i<NDIM;++i)
		invariant += dz[0][i]*dz[1][NDIM+i] - dz[1][i]*dz[0][NDIM+i];
	      invariant = fabs(invariant);
	      */
	 
	      if (invariant > score) 
		{
		  score = invariant;
		  kMax=k;		  
		  segment = sh[j];
		}	
	    }		
	}

      int segIndex=-1;
      // Find the segment that needs to be split by checking how much the highest
      // invariant will decrease depending on which segment is split

      //if ((computeSegment)&&(score>=threshold)) segIndex = length2<M>(nb.simplex,geometry).second; else 
      if ((computeSegment)&&(score>=threshold)) 
	{
	  int curSegIndex = nb.simplex->findSegmentIndex(segment);
	  double minScore = score;	  	  
	  Coord midPoint2[NDIM_W];	 

	  // First, try splitting the segment with maximal invariant value
	  // This gives two new invariants value, one for each midpoint of the split segment (=> coords[2])
	  coords[0]=nb.curV[kMax]->getCoordsPtr();
	  coords[1]=nb.neiV[kMax]->getCoordsPtr();	
	  coords[2]=midPoint2; // Midpoint of each part of the split segment
	  
	  // This is the midpoint of the maximal invariant segment
	  geometry->midPointCoords(segment->getVertex(0)->getCoordsPtr(),
				   segment->getVertex(1)->getCoordsPtr(),
				   midPoint);
	  
	  geometry->midPointCoords(midPoint,
				   segment->getVertex(1)->getCoordsPtr(),
				   midPoint2);
   
	  geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
	  double invariant1 = fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	 
	  geometry->midPointCoords(segment->getVertex(0)->getCoordsPtr(),
				   midPoint,
				   midPoint2);	  
   
	  geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
	  double invariant2 = fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	  
	  double invariant = std::max(invariant1,invariant2);
	  if (invariant<minScore)
	    {
	      minScore = invariant;
	      segIndex = curSegIndex;
	    }	  

	  // Next try splitting other segments => each split only gives one invariant, 
	  // computed by replacing the current simplex's vertex involved in the maximum 
	  //  invariant computation by the midPoint of the split segment

	  coords[0]=midPoint2; // The midpoint of the tentative segment (replaces current simplex's vertex)
	  coords[2]=midPoint;  // The midpoint corresponding to the maximaum value invariant (same as original)
	  for (int k=0; k<Simplex::NSEG; ++k)
	    {
	      if (k!=curSegIndex) // we already treated that particular case
		{		  
		  SegmentHandle sh = nb.simplex->getSegmentHandle(k);
		  // The segment must contain the vertex involved in the maximal invariant 
		  // computation (i.e. 'nb.curV[kMax]'). If it does not, then splitting it 
		  // will not affect the maximal invariant value !
		  if ((sh->getVertex(0)==nb.curV[kMax])||(sh->getVertex(1)==nb.curV[kMax]))
		    {
		      geometry->midPointCoords(sh->getVertex(0)->getCoordsPtr(),
					       sh->getVertex(1)->getCoordsPtr(),
					       midPoint2);
		      geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
		      invariant = fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
		      if (invariant<minScore)
			{
			  minScore = invariant;
			  segIndex = k;
			}	
		    }
		}
	    }
	  score -= minScore;
	}

      return CheckRefineReturnType(score,segIndex);
    }
        
    // -> Computes the (1-cosine) of the angle between the tracer, the vertices of 
    // a simplex and its barycenter.    
    template <class M>
    double
    tracerDisplacementAngle(const typename M::Simplex *simplex,
			    const typename M::GeometricProperties *geometry)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;
      typedef typename Simplex::Tracer::Type Tracer;

      //static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
     
      double score=1; // stores the cosine of an angle
      int maxIndex=0;
      //Simplex *simplex = nb.simplex;
      Tracer *tracer = simplex->tracer.getPointer(); 
      Coord barycenter[NDIM_W]={0};
      
      for (int i=0;i<Simplex::NVERT;++i) 
	{
	  const Coord *coords=simplex->getVertex(i)->getCoordsConstPtr();
	  
	  for (int j=0;j<NDIM_W;++j)
	    barycenter[j] +=geometry->
	      checkCoordConsistency(coords[j],tracer[j],j);
	}
      
      for (int i=0;i<NDIM_W;++i) 
	barycenter[i]/=Simplex::NVERT;
	
      for (int i=0;i<Simplex::NVERT;++i)
	{
	  Coord vec[2][NDIM_W];
	  geometry->getVector(simplex->getVertex(i)->getCoordsConstPtr(),
			      barycenter,vec[0]);			      
	  geometry->getVector(simplex->getVertex(i)->getCoordsConstPtr(),
			      tracer,vec[1]);

	  double n1=geometry->norm(vec[0]);
	  double n2=geometry->norm(vec[1]);
	  double tmpScore;
	  if ((n1>0)&&(n2>0))
	    {
	      tmpScore = geometry->dot_noCheck(vec[0],vec[1])/(n1*n2);
	    }
	  else tmpScore = 1;
	  
	  if (tmpScore<score) 
	    {
	      score=tmpScore;
	      maxIndex=i;
	    }
	}
      
      return 1.0L-score;      
    }

    template <class M>
    std::pair<double,int>
    poincareInvariantBelowSimplexTracer(typename M::Simplex* simplex,
					const typename M::GeometricProperties *geometry)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
  
      double score=0;
      int index=0;
      Coord midPoint[NDIM_W];
      Coord dz[2][NDIM_W];
      Coord *coords[3];

      /*
      Coord barycenter[NDIM_W]={0};      
      for (int i=0;i<Simplex::NVERT;++i) 
	{
	  const Coord *coords=nb.simplex->getVertex(i)->getCoordsConstPtr();
	  
	  for (int j=0;j<NDIM_W;++j)
	    barycenter[j] +=geometry->
	      checkCoordConsistency(coords[j],nb.simplex->tracer.getValueAt(j),j);
	}
      
      for (int i=0;i<NDIM_W;++i) 
	barycenter[i]/=Simplex::NVERT;
      */

      coords[1]=midPoint; // barycenter of opposite facet for each vertex k
      //coords[1]=barycenter;
      coords[2]=simplex->tracer.getPointer(); // tracer particle
      
      for (int k=0;k<Simplex::NVERT;++k)
	{
	  coords[0]=simplex->getVertex(k)->getCoordsPtr(); // vertex nb.curV[k]
	  std::fill_n(midPoint,NDIM_W,0);
	  
	  // Compute the barycenter of facet opposite to vertex k (->coords[1])
	  FacetHandle fh=simplex->getFacetHandle(k);
	  for (int j=0;j<Simplex::Facet::NVERT;++j)
	    {
	      const Coord* c = fh->getVertex(j)->getCoordsPtr(); //nb.fh[k]
	      for (int i=0;i<NDIM_W;++i)
		midPoint[i]+=geometry->checkCoordConsistency(c[i],coords[2][i],i);
	    }
	  
	  for (int i=0;i<NDIM_W;++i)
	    midPoint[i]/=Simplex::Facet::NVERT;
	  
	  // poincaré invariant computed from vertex/opp facet barycenter/tracer
	  geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
	  double invariant = fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));

	  if (invariant > score)
	    {
	      score = invariant;
	      index = k;
	    }
	}
      
      return std::make_pair(score,index);
    }

    template <class M>
    std::pair<double,int>
    poincareInvariantBelowSegTracers(typename M::Simplex* simplex,
				     //const typename M::SegmentHandle &seg, int sid,
				     const typename M::GeometricProperties *geometry,
				     bool average=false)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;
      
      Coord dz[2][NDIM_W];
      Coord *coords[3];
      double result=0;
      double avg=0;
      int index=0;

      for (int sid=0;sid<Simplex::NSEG;++sid)
	{
	  //if (sid!=excludeSeg) {
	    //int sid = seg->getSimplex()->findSegmentIndex(seg);
	    SegmentHandle seg=simplex->getSegmentHandle(sid);
	    coords[0]=seg->getVertex(0)->getCoordsPtr();
	    coords[1]=seg->getVertex(1)->getCoordsPtr();
	    coords[2]=&simplex->segTracers.getPointer()[sid*NDIM_W];
    
	    // poincaré invariant computed from vertex/opp facet barycenter/tracer
	    geometry->template getBaseVectors<Coord,Coord,2,NDIM_W>(coords,dz);
	    double pi=fabs(internal::computeInvariant<NDIM,NDIM_W>(dz));
	    if (pi>result) 
	      {
		result=pi;
		index=sid;
	      }
	    avg+=pi;
	    //}
	}
      if (average) return std::make_pair(avg/(Simplex::NSEG/*-(excludeSeg!=-1)*/),index);
      else return std::make_pair(result,index);
    }

    template <class M>
    std::pair<double,int>
    poincareInvariantBelowTracers(typename M::Simplex* simplex,
				  //const typename M::SegmentHandle &seg, int sid,
				  const typename M::GeometricProperties *geometry,
				  bool average=false)
    {
      typedef typename M::Simplex Simplex;

      std::pair<double,int> pi=
	poincareInvariantBelowSimplexTracer<M>(simplex,geometry);

      std::pair<double,int> spi=
	poincareInvariantBelowSegTracers<M>(simplex,geometry,average);

      if (pi.first>spi.first)
	{
	  if (average)
	    pi.first=
	      (spi.first*(Simplex::NSEG)+pi.first)/
	      (Simplex::NSEG+1);	      
	  
	  pi.second+=Simplex::NSEG;
	  return pi;	
	}
      else
	{
	  if (average)
	    spi.first=(spi.first*Simplex::NSEG+pi.first)/(Simplex::NSEG+1);	      
	  
	  return spi;	
	}
    }


    template <class M>
    typename M::CheckRefineReturnType 
    volumeBelow(const typename M::Simplex::Neighborhood &nb, 
		const typename M::GeometricProperties *geometry, 
		double threshold=0, bool computeSegment=true)
    {
      typedef typename M::CheckRefineReturnType CheckRefineReturnType;
      typedef typename M::SegmentHandle SegmentHandle;
      typedef typename M::FacetHandle   FacetHandle;
      typedef typename M::Simplex Simplex;   
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord   Coord;

      static const int NDIM = M::NDIM;
      static const int NDIM_W = M::NDIM_W;

      if (nb.nNeighbors!=Simplex::NNEI) return CheckRefineReturnType(0,-1);
    
      Coord *vCoordsPtr[NDIM_W+1];
      Coord base[NDIM_W][NDIM_W];
      double volume[NDIM+1];

      for (int i=0;i<Simplex::NVERT;++i) 
	vCoordsPtr[i]=nb.curV[i]->getCoordsPtr();
	 
      // Check the volume of each of the (NDIM-1) NDIM_W-simplices below 's' 
      // that we can build with the (NDIM+1) vertices 'curV 'of the current simplex and
      // the (NDIM+1) vertices 'neiV' that belong to neighboring simplices
      // Each NDIM_W-simplices is made of the vertices of 's' and opposite vertices
      // of (NDIM) of the (NDIM+1) vertices of the neighbors.
      for (int j=Simplex::NVERT;j<=NDIM_W+1;++j)
	{		
	  int k=0;
	  for (int i=Simplex::NVERT;i<NDIM_W+1;++i,++k) 
	    {
	      if (i==j) ++k;
	      vCoordsPtr[i]=nb.neiV[k]->getCoordsPtr();
	    }

	  geometry->template getBaseVectors<Coord,Coord,NDIM_W,NDIM_W>(vCoordsPtr,base);

	  // The (j-Simplex::NVERT)th neighbor is not involved in 'volume[j-Simplex::NVERT]'
	  volume[j-Simplex::NVERT] = 
	    SimplexVolumeT<NDIM_W,NDIM_W,Coord>::compute(base);
	}
	
      int maxId=0;
      
      for (int i=1;i<NDIM+1;++i)
	if (volume[maxId]<volume[i]) maxId=i;
      
      double score=volume[maxId];
      int segIndex = -1;
      if ((computeSegment)&&(score>=threshold)) 
	segIndex = length2<M>(nb.simplex,geometry).second;
	
      return CheckRefineReturnType(score,segIndex);
    }

    template <class M>
    class CellDataFunctor_lengthT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_lengthT(const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	return sqrt(length2<M>(const_cast<Cell*>(c),Base::mesh->getGeometry(),0,false).first);
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("longestSide");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_poincareInvariantT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_poincareInvariantT(const M* m,int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	//typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return poincareInvariant<M>
	  (const_cast<Cell*>(c),Base::mesh->getGeometry()).first;
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("poincareInvariant");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_poincareInvariantWithSegTracersT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_poincareInvariantWithSegTracersT
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	//typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return poincareInvariantWithSegTracers<M>
	  (const_cast<Cell*>(c),Base::mesh->getGeometry());
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("poincareInvariantWithSegTracers");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_poincareInvariantWithSegTracers_order1T: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_poincareInvariantWithSegTracers_order1T
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	//typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return poincareInvariantWithSegTracers_order1<M>
	  (const_cast<Cell*>(c),Base::mesh->getGeometry()).first;
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const 
      {return std::string("poincareInvariantWithSegTracers_O1");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_poincareInvariantFromNeighborsT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_poincareInvariantFromNeighborsT
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return poincareInvariantFromNeighbors<M>
	  (nb,Base::mesh->getGeometry(),0,false).first;
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("poincareInvariantFromNeighbors");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_poincareInvariantBelowTracersT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_poincareInvariantBelowTracersT
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	//typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return poincareInvariantBelowTracers<M>
	  (const_cast<Cell*>(c),Base::mesh->getGeometry()).first;
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("poincareInvariantBelowTracers");}
      int getSize() const {return 1;}
    }; 
    
    template <class M>
    class CellDataFunctor_tracerDisplacementAngleT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_tracerDisplacementAngleT
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	//typename Cell::Neighborhood nb( const_cast<Cell*>(c) ); 
	return tracerDisplacementAngle<M>
	  (const_cast<Cell*>(c),Base::mesh->getGeometry());
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("tracerDisplacementAngle");}
      int getSize() const {return 1;}
    }; 

    template <class M>
    class CellDataFunctor_volumeBelowT: 
      public cellDataFunctors::InterfaceT<M,typename M::Simplex>
    {
    public:
      typedef cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
      typedef typename Base::Cell Cell;         
      CellDataFunctor_volumeBelowT
      (const M* m, int flags=cellDataFunctors::F_NO_FLAG):
	Base(m,flags){}

      double get(const Cell *c) const
      {
	// FIXME: const_cast => Design flaw
	typename Cell::Neighborhood nb( const_cast<Cell*>(c) );
	return volumeBelow<M>(nb,Base::mesh->getGeometry(),0,false).first;
      }

      double get(const Cell *c, int n) const {return get(c);}
      std::string getName() const {return std::string("volumeBelow");}
      int getSize() const {return 1;}
    }; 
    
  } // namespace refine
} // namespace slv

/** \}*/

#include "../internal/namespace.footer"
#endif
