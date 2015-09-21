#ifndef __LOCAL_AMR_GRID_PROJECTOR_BASE_2D__HXX__
#define __LOCAL_AMR_GRID_PROJECTOR_BASE_2D__HXX__

#include "../../tools/wrappers/boostMultiprecisionFloat128.hxx"
#include "localAmrGridProjectorBasePrototype.hxx"

#include "../../internal/namespace.header"

namespace internal {

  namespace {
    template <int NDIM, class M, typename F, typename HF>
    struct FaceVectors;
    
    // Generic version
    template <class M, typename F>
    struct FaceVectors<2,M,F,F>
    {
      static const int NDIM=2;

      typedef F Float;
      typedef F HFloat;      

      typedef typename M::Simplex Simplex;
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord Coord;

      template <class Edge, class G>
      static void computeT(Edge &edge, G *meshGeometry, Float tolerance=1.E-3)
      {
	meshGeometry->template
	  getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
				      edge.otherVertex->getCoordsConstPtr(),
				      edge.vecT);
		  	
	Float norm=meshGeometry->template normalize_noCheck<Float,NDIM>(edge.vecT);
	if (norm==0) 
	  {		  
	    PRINT_SRC_INFO(LOG_WARNING);
	    glb::console->print<LOG_WARNING>("I encountered a highly degenerate triangle (two vertices at the exact same location). Cross your fingers ;)\n");
	    edge.vecT[0]=1;edge.vecT[1]=0;
	  }
      }      
    };

     // This version uses HF only if F is not precise enough
    template <class M, typename F, typename HF>
    struct FaceVectors<2,M,F,HF>
    {
      static const int NDIM=2;

      typedef F Float;
      typedef HF HFloat;      

      typedef typename M::Simplex Simplex;
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord Coord;

      template <class Edge, class G>
      static void computeT(Edge &edge, G *meshGeometry, Float tolerance=1.E-3)
      {
	FaceVectors<NDIM,M,F,F>::computeT(edge,meshGeometry,tolerance);

	if ((fabs(edge.vecT[0])<tolerance)||
	    (fabs(edge.vecT[1])<tolerance))
	  {
	    HFloat vecHT[NDIM];
	    meshGeometry->template
	      getVector<Coord,HFloat,NDIM>(edge.vertex->getCoordsConstPtr(),
					   edge.otherVertex->getCoordsConstPtr(),
					   vecHT);
		  
	    HFloat norm=meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHT);
		 
	    if (norm==0) 
	      {		  
		PRINT_SRC_INFO(LOG_WARNING);
		glb::console->print<LOG_WARNING>("I encountered a highly degenerate triangle (two vertices at the exact same location). Cross your fingers ;)\n");
		edge.vecT[0]=1;
		edge.vecT[1]=0;
	      }
	    else
	      {
		edge.vecT[0]=hlp::numericStaticCast<Float>(vecHT[0]);
		edge.vecT[1]=hlp::numericStaticCast<Float>(vecHT[1]);		    
	      }
	  }
      }
    };
  }

  template <class AMR, class MESH, typename F, typename HF, bool checkAccuracy, int dummy> 
  class LocalAmrGridProjectorBaseT<2,AMR,MESH,F,HF,checkAccuracy,dummy>
  {
  public:
    static const int NDIM = 2;

    typedef typename MESH::vertexPtr_iterator vertexPtr_iterator;
    typedef typename MESH::simplexPtr_iterator simplexPtr_iterator;
    typedef typename MESH::Simplex Simplex;
    typedef typename MESH::Vertex  Vertex;
    typedef typename MESH::Segment Segment;
    typedef typename MESH::SegmentHandle SegmentHandle;
    typedef typename MESH::Coord Coord;
    typedef typename MESH::GeometricProperties MeshGeometricProperties;

    typedef typename AMR::ICoord ICoord;
    typedef typename AMR::Voxel Voxel;
    typedef typename AMR::GeometricProperties AmrGeometricProperties;

    typedef F Float;
    typedef HF HFloat;
    // FILE *fl;

    LocalAmrGridProjectorBaseT(AMR *amr_, MESH *mesh_, 
			       int nThreads_)://,MpiCommunication *mpiCom_):
      amr(amr_),
      mesh(mesh_),  
      nThreads(nThreads_)//,
      //mpiCom(mpiCom_)
    {
      meshGeometry = mesh->getGeometry();
      amrGeometry = amr->getGeometry();
      // fl=fopen("vert.dat","w");      
    }

    virtual ~LocalAmrGridProjectorBaseT()
    {
      //fclose(fl);
    }

 
  protected:
    AMR *amr;
    MESH *mesh;  
    int nThreads; 
    //MpiCommunication *mpiCom;

    MeshGeometricProperties *meshGeometry;
    AmrGeometricProperties *amrGeometry;
    
  public:
    struct IncidentEdge {
      static const int flag_boundary = (1<<0);
      static const int flag_fold = (1<<1);

      Simplex *simplex;       // current simplex
      Vertex  *vertex;        // current vertex 

      Simplex *otherSimplex;  // other simplex (i.e. opposite to simplex wrt edge)
      Vertex  *otherVertex;   // other vertex (i.e. other extremity of the segment)      

      //Float  weight;       // weight of the current simplex	
      //Float  otherWeight;  // weight of the other simplex
      //Float foldingFactor; // {+1 or -1} when the edge is normal / a fold  
      //Float otherContrib;  // Contrib to the other extremity of the edge

      Float deltaWeight;
      Float deltaGrad[NDIM];      
      Float vecT[NDIM];        // Tangent vector
      Float vecN[NDIM];	       // Normal vector;
      int flags;

      void print()
      {
	printf ("S:%ld V:%ld OS:%ld OV:%ld\n",simplex,vertex,otherSimplex,otherVertex);
	printf("deltaW = %e\n",(double)deltaWeight);
	printf("deltaGrad = %e\n",(double)deltaGrad[0],(double)deltaGrad[1]);
	printf("vecT = %e\n",(double)vecT[0],(double)vecT[1]);
	printf("vecN = %e\n",(double)vecN[0],(double)vecN[1]);
	printf("flags = %d\n",flags);
	simplex->template print<LOG_STD>("S:");
	otherSimplex->template print<LOG_STD>("OS:");
	vertex->template print<LOG_STD>("");
	otherVertex->template print<LOG_STD>("");
      }
    };

    template <class WF>
    void getVertexContribAndEdges(Simplex * const * simplex, int nSimplices,
				  Vertex *vertex, const WF &wf,
				  std::vector<IncidentEdge> &ownedEdges,
				  bool debug=false)
    {
      ownedEdges.clear();
  
      for (int i=0;i<nSimplices;++i)
	getEdges(simplex[i],vertex,wf,ownedEdges);
      /*
      Float contrib=0;
      for (unsigned long i=0;i<ownedEdges.size();++i)
	contrib += computeContrib(vertex->getCoordsConstPtr(),wf,ownedEdges[i]);
      
      return contrib;
      */
    }

    // FIXME : Should we always compute contributions in the coordinates of the voxel
    // to avoid rounding errors when coordinates are large ?
    template <class WF>
    void getEdges(Simplex *simplex, Vertex *vertex, const WF& wf,
		  std::vector<IncidentEdge> &ownedEdges)
		  
    {
      bool checked[2];
      IncidentEdge e[2];
      int          owned[2];
      Float        tmpVec[NDIM];
      Float        weight;
      Float        grad[NDIM];
      int          nOwned=0;

      if (invalidSimplex(simplex)) return;

      //weight=simplex->cache.d;
      char tag=wf.getTag(simplex);
      weight = wf.template get<Float>(simplex);
      if (WF::ORDER>0) wf.getGradient(simplex,grad);

      // Get the two segments incident to 'vertex' as well as the corresponding
      // neighboring simplex
      if (simplex->getVertex(0)==vertex) 
	{
	  e[0].otherVertex=simplex->getVertex(1);
	  e[0].otherSimplex=simplex->getNeighbor(2);
	  e[1].otherVertex=simplex->getVertex(2);
	  e[1].otherSimplex=simplex->getNeighbor(1);		
	}
      else if (simplex->getVertex(1)==vertex) 
	{
	  e[0].otherVertex=simplex->getVertex(0);
	  e[0].otherSimplex=simplex->getNeighbor(2);
	  e[1].otherVertex=simplex->getVertex(2);
	  e[1].otherSimplex=simplex->getNeighbor(0);
	}
      else
	{
	  e[0].otherVertex=simplex->getVertex(0);
	  e[0].otherSimplex=simplex->getNeighbor(1);
	  e[1].otherVertex=simplex->getVertex(1);
	  e[1].otherSimplex=simplex->getNeighbor(0);
	}
     
      if (tag)
	{
	  checked[0]=
	    (invalidSimplex(e[0].otherSimplex))||
	    ((e[0].otherSimplex<simplex)&&(!wf.getTag(e[0].otherSimplex)));

	  checked[1]=
	    (invalidSimplex(e[1].otherSimplex))||
	    ((e[1].otherSimplex<simplex)&&(!wf.getTag(e[1].otherSimplex)));
	  /*
	  otherTag[0]=invalidSimplex(e[0].otherSimplex)?0:wf.getTag(e[0].otherSimplex);
	  otherTag[1]=invalidSimplex(e[1].otherSimplex)?0:wf.getTag(e[1].otherSimplex);
	  */
	}
      else
	{	  
	  checked[0]=
	    (invalidSimplex(e[0].otherSimplex))||
	    (e[0].otherSimplex<simplex);

	  checked[1]=
	    (invalidSimplex(e[1].otherSimplex))||
	    (e[1].otherSimplex<simplex);	 
	}

      // We do not want to compute the contribution of each segment twice
      // so we only consider a segment owned if 'otherSimplex' has an address
      // lower than 'simplex' (or else it will be computed when 'otherSimplex'
      // is the current simplex). Note that it still works when otherSimplex
      // is NULL (i.e. on the boundary of a mesh)
      // FIXME : Check if the volume of the simplex is non 0 to see 
      // if it is degenerate ?
      for (int j=0;j<2;j++)
	{
	  /*
	  if ((e[j].otherVertex<vertex)&&
	      ((e[j].otherSimplex<simplex)||invalidSimplex(e[j].otherSimplex)))	  
	  */
	  if ((e[j].otherVertex<vertex)&&(checked[j]))
	    {	
	      // Complete missing data
	      e[j].simplex=simplex;
	      e[j].vertex=vertex;

	      // Compute tangent and normal vectors
	      FaceVectors<NDIM,MESH,Float,HFloat>::computeT(e[j],meshGeometry);
	      /* // this was replaced by computeT above	      
	      meshGeometry->template
		getVector<Coord,Float,NDIM>(vertex->getCoordsConstPtr(),
					    e[j].otherVertex->getCoordsConstPtr(),
					    e[j].vecT);
		  
	      // FIXME: check if norm==0 ?	      
	      Float tNorm=meshGeometry->template normalize_noCheck<Float,NDIM>(e[j].vecT);
	      if (tNorm==0) 
		{		  
		  PRINT_SRC_INFO(LOG_WARNING);
		  glb::console->print<LOG_WARNING>("I encountered a highly degenerate triangle (two vertices at the exact same location). Cross your fingers ;)\n");
		  e[j].vecT[0]=1;e[j].vecT[1]=0;
		}

	      if ((fabs(e[j].vecT[0])<1.E-3)||
		  (fabs(e[j].vecT[1])<1.E-3))
		{
		  HFloat vecHT[NDIM];
		  meshGeometry->template
		    getVector<Coord,HFloat,NDIM>(vertex->getCoordsConstPtr(),
						 e[j].otherVertex->getCoordsConstPtr(),
						 vecHT);
		  
		  HFloat normH=meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHT);
		 
		  if (normH==0) 
		    {		  
		      PRINT_SRC_INFO(LOG_WARNING);
		      glb::console->print<LOG_WARNING>("I encountered a highly degenerate triangle (two vertices at the exact same location). Cross your fingers ;)\n");
		      e[j].vecT[0]=1;
		      e[j].vecT[1]=0;
		    }
		  else
		    {
		      e[j].vecT[0]=hlp::numericStaticCast<Float>(vecHT[0]);
		      e[j].vecT[1]=hlp::numericStaticCast<Float>(vecHT[1]);		    
		    }
		}
	      */

	      // Normal vector is orthogonal to tangent vector
	      e[j].vecN[0] = -e[j].vecT[1];
	      e[j].vecN[1] =  e[j].vecT[0];	      	      
	      
	      /*
	      e[j].weight=weight;

	      if (validSimplex(e[j].otherSimplex))
		e[j].otherWeight=e[j].otherSimplex->cache.d;
	      else
		e[j].otherWeight=0;
	      */

	      e[j].flags=0;

	      if (validSimplex(e[j].otherSimplex))
		{
		  e[j].deltaWeight=wf.template get<Float>(e[j].otherSimplex);
		  if (WF::ORDER>0) wf.getGradient(e[j].otherSimplex,e[j].deltaGrad);
		}
	      else
		{
		  e[j].deltaWeight=0;
		  e[j].flags|=IncidentEdge::flag_boundary; // Boundary edge

		  if (WF::ORDER>0)
		    {
		      e[j].deltaGrad[0]=0;
		      e[j].deltaGrad[1]=0;
		    }
		}
	      

	      /*
	      if ((e[j].otherSimplex==NULL)||(e[j].otherSimplex->isGhost())) 
		e[j].otherWeight=0;
	      else 
		e[j].otherWeight=e[j].otherSimplex->cache.d;
	      */

	      //tmpVec is the other segment in this simplex departing from vertex
	      meshGeometry->template
		getVector<Coord,Float,NDIM>(vertex->getCoordsConstPtr(),
					    e[(j+1)%2].otherVertex->getCoordsConstPtr(),
					    tmpVec);
	      
	      // Check if the normal points inside or outside the current simplex
	      // and negate it if it points outside
	      // FIXME: check if dot==0 ? In principle that should be fine even 
	      // for degenerate cases ...

	      //meshGeometry->template normalize_noCheck<Float,NDIM>(tmpVec);
	      Float dot = meshGeometry->template
		dot_noCheck<Float,Float,NDIM>(e[j].vecN,tmpVec);
	      //if (fabs(dot)<1.E-16) {fprintf(stderr,"DOT=0\n");exit(-1);} else
	      if (dot<0)
		{
		  e[j].vecN[0]=-e[j].vecN[0];
		  e[j].vecN[1]=-e[j].vecN[1];		  
		}	      

	      // Test whether the sheet is folding along this edge 
	      if (e[j].otherSimplex!=NULL)
		{
		  /*
		  int neiId = simplex->getNeighborIndex(e[j].otherSimplex);
		  e[j].foldingFactor = 
		    simplex->compareNeighborOrientation(neiId);
		  */
		  
		  int vid   = e[j].otherSimplex->getNeighborIndex(simplex);
		  Vertex *v = e[j].otherSimplex->getVertex(vid);
		  		 
		  meshGeometry->template
		    getVector<Coord,Float,NDIM>(vertex->getCoordsConstPtr(),
						v->getCoordsConstPtr(),
						tmpVec);

		  //meshGeometry->template normalize_noCheck<Float,NDIM>(tmpVec);
		  Float dot2 = meshGeometry->template
		    dot_noCheck<Float,Float,NDIM>(e[j].vecN,tmpVec);
		  
		  // Don't forget vecN has been negated if (dot<0)

		   // if (dot2>0) e[j].foldingFactor=-1;
		   // else e[j].foldingFactor=1;
		   
		  
		  // The sheet is folding here, so the neighbor
		  // triangles have their common face normal pointing
		  // in the same direction
		  
		  //if (fabs(dot2)<1.E-16) {fprintf(stderr,"DOT2=0\n");exit(-1);} else
		  if (dot2>0) 
		    {
		      e[j].deltaWeight = -e[j].deltaWeight;		  
		      if (WF::ORDER>0)
			{
			  e[j].deltaGrad[0] = -e[j].deltaGrad[0];
			  e[j].deltaGrad[1] = -e[j].deltaGrad[1];
			  e[j].flags|=IncidentEdge::flag_fold; // fold edge
			}
		    }
       
		}	      

	      e[j].deltaWeight = weight - e[j].deltaWeight;
	      if (WF::ORDER>0)
		{		  
		  e[j].deltaGrad[0] = grad[0] - e[j].deltaGrad[0];
		  e[j].deltaGrad[1] = grad[1] - e[j].deltaGrad[1];
		}
	      owned[nOwned++]=j; //Mark as owned by the simplex and vertex
	    }
	}

      for (int j=0;j<nOwned;++j)
	ownedEdges.push_back(e[owned[j]]);
      /*
      Float contrib=0;
      // Now compute the vertex contribution if needed
      // And push_back the edges owned by the current vertex AND simplex
      if (nOwned)
	{	  
	  for (int j=0;j<nOwned;++j)
	    {
	      IncidentEdge &edge=e[owned[j]];
	      contrib += computeContrib(vertex->getCoordsConstPtr(),edge);
	      //edge.otherContrib = -computeContrib(edge.otherVertex->getCoordsConstPtr(),edge);	    
	      ownedEdges.push_back(edge);
	    }	   
	}
    
      return contrib;
      */
    }

    // Contribution to the current voxel due to mesh segments intersecting voxel faces .    
    template <typename R, class WF, class OutputIterator>
    int addVoxelFacetSimplexEdgeContrib(const R &raytracer, const WF &wf, 
					const IncidentEdge &edge, OutputIterator out)
    {
    
      typedef typename R::Coord RCoord;
      const RCoord *exitPoint=raytracer.getExitPointConstPtr();

      int normalDim = raytracer.getExitFaceNormalDim();     
   
      // Sign of the Tangent vector in the exitPoint/voxelFace contrib
      Float Tsign = hlp::sgn(edge.vecN[1-normalDim]);
      // Sign of the Normal vector in the exitPoint/voxelFace contrib
      Float Nsign = raytracer.getExitFaceNormalSign();      
      
      Voxel *curVoxel = raytracer.getCurVoxel();
      Voxel *nextVoxel = raytracer.getNextVoxel();
      Float contrib;
      
      
      // This happens if the simplex edge is tengant to the voxel edge it crosses
      //bool almostDegenerate=(fabs(edge.vecT[normalDim])<1.E-10);      
      if (Tsign==0)
      //if (fabs(edge.vecT[normalDim])<1.E-10)
	{
	  // This is needed because vecT may not be precise enough to determine the size
	  // as it is already a difference between coordinates !
	  Tsign = meshGeometry->coordDiffSign(edge.vertex->getCoord(normalDim),
					      edge.otherVertex->getCoord(normalDim),
					      normalDim);
	  //printf("Tsign = %f\n",(double)Tsign);
	  if ((Tsign>0)==(edge.vecN[normalDim]>0))
	    Tsign=-hlp::sgn(edge.vecT[1-normalDim]);
	  else
	    Tsign=hlp::sgn(edge.vecT[1-normalDim]);
	}

      if (AMR::BOUNDARY_TYPE != BoundaryType::PERIODIC)
	{
	  /*
	  // There are 2 contrib from the exitPoint / voxel face    
	  contrib=computeOrthogonalContrib(exitPoint,Tsign,Nsign,normalDim,wf,edge);
	  // Plus 2 contribs from the exitPoint / ray
	  contrib+=computeContrib(exitPoint,wf,edge);
	  */
	  contrib=computeVoxelFacetSimplexEdgeContrib(exitPoint,Tsign,Nsign,normalDim,
						      wf,edge);

	  // Contribution to each voxel is opposite
	  //contrib-=maxContrib;
	  (*out)=std::make_pair(nextVoxel,contrib);++out;
	  (*out)=std::make_pair(curVoxel,-contrib);++out;
	  // (*out)=std::make_pair(nextVoxel,maxContrib);++out;
	  // (*out)=std::make_pair(curVoxel,-maxContrib);++out;
	}
      else
	{
	  // if coordinates are on a periodic boundary, we need to shift them to the 
	  // same side as the intersection polygon they will contribute to !	 
	  int boundary=0;

	  for (int i=0;i<NDIM;++i) 
	    {
	      if (amr->getBBoxMax(i)-amr->getEpsilon(i)<exitPoint[i]) 
		boundary|=1<<i;
	      else if (amr->getBBoxMin(i)+amr->getEpsilon(i)>exitPoint[i]) 
		boundary|=1<<i;
	    }
	 
	  // if (check) printf("OnBoundary = %d\n",boundary);
	  if (boundary)
	    {
	      // Coordinates are on a boundary, so their contribution on either side will 
	      // be different as we have to shift them so that they are close to the voxel
	      // they will contribute to !
	      //printf("BOUNDARY!\n");
	      RCoord refPoint[NDIM];
	      //RCoord refPoint_Next[NDIM];
	      RCoord exitPoint_Cur[NDIM];
	      RCoord exitPoint_Next[NDIM];

	      amr->index2CenterCoords(curVoxel->getIndex(),refPoint);
	      for (int i=0;i<NDIM;++i)
		exitPoint_Cur[i]=refPoint[i]+ amrGeometry->template 
		  correctCoordsDiff<RCoord>(exitPoint[i]-refPoint[i],i);

	      amr->index2CenterCoords(nextVoxel->getIndex(),refPoint);
	      for (int i=0;i<NDIM;++i)
		exitPoint_Next[i]=refPoint[i]+ amrGeometry->template 
		  correctCoordsDiff<RCoord>(exitPoint[i]-refPoint[i],i);

	      
	      // if (raytracer.getCurVoxel()->getIndex()==4602678820611293184L)
	      // 	printf("Cur =[%lg,%lg], next = [%lg,%lg], Tsign =%f\n",
	      // 	       (double)exitPoint_Cur[0],(double)exitPoint_Cur[1],
	      // 	       (double)exitPoint_Next[0],(double)exitPoint_Next[1],(double)Tsign);

	      // if (check) printf("Cur =[%g,%g], next = [%g,%g] %g\n",
	      // 			exitPoint_Cur[0],exitPoint_Cur[1],
	      // 			exitPoint_Next[0],exitPoint_Next[1],
	      // 			edge.foldingFactor);

	      // Contrib to each voxel is compute separately
	      /*
	      Float contribCur=computeContrib(exitPoint_Cur,wf,edge);
	      contribCur+=
		computeOrthogonalContrib(exitPoint_Cur,Tsign,Nsign,normalDim,wf,edge);
	      */
	      //std::cout<<"CUR: ("<<exitPoint_Cur[0]<<","<<exitPoint_Cur[1]<<")\n";
	      Float contribCur=
		computeVoxelFacetSimplexEdgeContrib(exitPoint_Cur,Tsign,Nsign,normalDim,
						    wf,edge);
	      //contribCur-=maxContrib;
	      (*out)=std::make_pair(curVoxel,-contribCur);++out;
	      //(*out)=std::make_pair(curVoxel,-maxContrib);++out;

	      /*
	      Float contribNext=computeContrib(exitPoint_Next,wf,edge);
	      contribNext+=
		computeOrthogonalContrib(exitPoint_Next,Tsign,Nsign,normalDim,wf,edge);
	      */
	      //std::cout<<"NXT: ("<<exitPoint_Next[0]<<","<<exitPoint_Next[1]<<")\n";
	      Float contribNext=
		computeVoxelFacetSimplexEdgeContrib(exitPoint_Next,Tsign,Nsign,normalDim,
						    wf,edge);
	      //contribNext-=maxContrib;
	      (*out)=std::make_pair(nextVoxel,contribNext);++out;
	      //(*out)=std::make_pair(nextVoxel,maxContrib);++out;
	      
	      contrib = contribCur;
	    }
	  else
	    {
	      /*
	      contrib=computeContrib(exitPoint,wf,edge);
	      contrib+=computeOrthogonalContrib(exitPoint,Tsign,Nsign,normalDim,wf,edge);
	      */
	      contrib =
		computeVoxelFacetSimplexEdgeContrib(exitPoint,Tsign,Nsign,normalDim,
						    wf,edge);

	      // Contribution to each voxel is opposite
	      //contrib-=maxContrib;
	      (*out)=std::make_pair(nextVoxel,contrib);++out;
	      (*out)=std::make_pair(curVoxel,-contrib);++out;
	      //(*out)=std::make_pair(nextVoxel,maxContrib);++out;
	      //(*out)=std::make_pair(curVoxel,-maxContrib);++out;
	    }
	  // fprintf(fl,"%g %g %g\n",exitPoint[0],exitPoint[1],contrib);
	}

      return 2;     
    }

    template <class WF, class OutputIterator>
    int addVoxelEdgeSimplexFacetContrib(Simplex *simplex, 
					int simplexNeighborIndex, 
					Float projectedFacetNormal[NDIM], int dim,
					Float contribSign[NDIM], 
					Float intersection[NDIM],
					Voxel *segNeighbor[Voxel::NSEGNEI],	
					const WF &wf, OutputIterator out,
					bool coordsAreConsistent=false,
					bool dbg=false) const    
    {
      return 0;
    }  
    
  protected:

    bool validSimplex(Simplex *simplex) const
    {
      if ((simplex==NULL)||(!simplex->isLocal()))
	return false;
      else 
	return true;
    }

    bool invalidSimplex(Simplex *simplex) const
    {
      if ((simplex==NULL)||(!simplex->isLocal()))
	return true;
      else 
	return false;
    }

    template <int sign, class CT, class WF, class OutputIterator>
    int addContrib(const CT  * coord, const WF &wf, const IncidentEdge &edge,
		   Voxel *v, OutputIterator out, bool debug=false)
    {
      Float contrib=computeContrib(coord,wf,edge,debug);

      if (sign>0)
	(*out)=std::make_pair(v,contrib);
      else
	(*out)=std::make_pair(v,-contrib);

      ++out;	  

      return 1;
    }

    // Compute the 2 contribs from a point ray intersection, where the tangent vectors of
    // adjacent simplices are the same, but their normal vectors are opposite
    // Contribution of each segment is w_i * P.T_i * P.N_i (order 0)
    // w_i weight of the simplex
    // P Coordinate of the vertex
    // T_i normalized vector tangeant to the segment
    // N_i normalized vector normal to the segment    
    template <class T, class WF>
    Float computeContrib(const T  * coord, const WF &wf, const IncidentEdge &edge,
			 bool debug=false)
    {         
      Float contrib = edge.deltaWeight;
      Float PT = edge.vecT[0]*coord[0] + edge.vecT[1]*coord[1];
      Float PN = edge.vecN[0]*coord[0] + edge.vecN[1]*coord[1];
      
      if (WF::ORDER>0)
	{
	  Float GT=edge.vecT[0]*edge.deltaGrad[0]+edge.vecT[1]*edge.deltaGrad[1];
	  Float GN=edge.vecN[0]*edge.deltaGrad[0]+edge.vecN[1]*edge.deltaGrad[1];
	  
	  if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
	    contrib = getPeriodizedDeltaWeight(coord,wf,edge);
	 	  
	  contrib += (GN*PN*2.0 + GT*PT)/3.0;	  
	}
      
      return contrib*PT*PN;
    }

    /*
    // Compute the two contributions (one for each voxel) from an axis aligned plane
    // (voxel facet) intersection with a ray (simplex edge). The tangent vectors of
    // each contribution is opposite, while their normal contrib is the same.
    // TSign and NSign are the signs (i.e. +1 or -1) of the Tangent and Normals.
    template <class T, class T2, class WF>
    Float computeOrthogonalContrib(const T *coord, T2 Tsign, T2 Nsign, int normalDim,
				    const WF &wf, const IncidentEdge &edge)
    {
      Float contrib = edge.deltaWeight;
      Float PTPN = coord[0]*coord[1]*Tsign*Nsign;
    
      if (WF::ORDER>0)
	{	  
	  Float GNPN2=edge.deltaGrad[normalDim]*coord[normalDim]*2.0; 
	  Float GTPT=edge.deltaGrad[1-normalDim]*coord[1-normalDim]; 

	  if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
	    contrib = getPeriodizedDeltaWeight(coord,wf,edge);

	  contrib += (GTPT + GNPN2)/3.0;	 
	}
      
      return contrib*PTPN;
    }
    */

    template <class T, class T2, class WF>
    Float computeVoxelFacetSimplexEdgeContrib(const T  * coord, 
					      T2 Tsign, T2 Nsign, int normalDim,
					      const WF &wf, const IncidentEdge &edge)
    {
      Float contrib1 = edge.deltaWeight; // Contrib along simplex edge
      Float contrib2 = edge.deltaWeight; // Contrib orthogonal to voxel face
      Float PT1 = edge.vecT[0]*coord[0] + edge.vecT[1]*coord[1];
      Float PN1 = edge.vecN[0]*coord[0] + edge.vecN[1]*coord[1];
      Float PTPN2 = coord[0]*coord[1]*Tsign*Nsign;

      // std::cout.precision(20);
      // std::cout <<std::scientific<<contrib1 << " " <<contrib2 <<" " <<PT1 <<" " 
      // 		<<PN1 << " " <<PTPN2 <<std::endl;
      
      if (WF::ORDER>0)
	{
	  Float GT1=edge.vecT[0]*edge.deltaGrad[0]+edge.vecT[1]*edge.deltaGrad[1];
	  Float GN1=edge.vecN[0]*edge.deltaGrad[0]+edge.vecN[1]*edge.deltaGrad[1];
	  Float GNPN2=edge.deltaGrad[normalDim]*coord[normalDim]; 
	  Float GTPT2=edge.deltaGrad[1-normalDim]*coord[1-normalDim]; 

	  if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
	    contrib1=contrib2=getPeriodizedDeltaWeight(coord,wf,edge);
	 
	  contrib1 += (GN1*PN1*2.0 + GT1*PT1)/3.0;	
	  contrib2 += (GNPN2*2.0 + GTPT2)/3.0;	 
	}
      
      return contrib1*PT1*PN1 + contrib2*PTPN2;
    }

    template <class T, class T2, class WF,class OutputIterator>
    int addVoxelCornerContrib(T *coords, 
			      T2 vertexNeighborDir[1<<NDIM][NDIM],
			      Voxel *neighbors[1<<NDIM],
			      const WF &wf, 
			      Simplex *simplex,
			      OutputIterator out, 
			      bool coordsAreConsistent=false)
    {
      // We add the contribution to each incident voxel.
      // Note that some incident voxels should not get a contribution,
      // as they do not share the vertex (it belongs to their face),
      // but we can still count them, as their contributions will 
      // cancel out.
      static const Float dimFactorInv = 2.0;
      
      Float contrib0=dimFactorInv*
	getPeriodizedWeight(coords,wf,simplex);//,coordsAreConsistent);	
      Float grad[NDIM];
      if (WF::ORDER>0) wf.getGradient(simplex,grad);

      T newCoords[NDIM];

      for (unsigned long n=0;n<(1<<NDIM);++n) 
	{
	  Float contrib=0;
	  Float PTPN=1.0;
	  int noShiftCount=0;
	  
	  for (int d=0;d<NDIM;++d) // each dim
	    {		
	      T c=coords[d];
	      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
		{
		  // if coordinates are on a periodic boundary, we need
		  // to shift them to the same side as the intersection
		  // polygon they will contribute to !
		  if ((c>=amr->getBBoxMax(d)-amr->getEpsilon(d))&&
		      (vertexNeighborDir[n][d]>0))
		    c=amr->getBBoxMin(d);
		  else if ((c<=amr->getBBoxMin(d)+amr->getEpsilon(d))&&
			   (vertexNeighborDir[n][d]<0))
		    c=amr->getBBoxMax(d);
		  else noShiftCount++;

		  newCoords[d]=c;
		}

	      Float PV=c*vertexNeighborDir[n][d];
	      if (WF::ORDER>0)
		{
		  Float GV=vertexNeighborDir[n][d]*grad[d];
		  contrib += PV*GV;
		}
	      PTPN *= PV;	    
	    }
			   
	  if (WF::ORDER>0)
	    {
	      // For periodic boundary crossing simplices, the gradient
	      // contrib may need update if the coordinates of a voxel 
	      // corner lying on the edge is shifted to any equivalent 
	      // coordinates on the boundary (i.e. if noShiftCount!=NDIM)
	      
	      if ((AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)&&
		  (noShiftCount!=NDIM))
		{
		  contrib+= dimFactorInv*
		    getPeriodizedWeight(newCoords,wf,simplex)-
		    contrib0;		  
		}
				
	      contrib=PTPN*(contrib0+contrib);
	    }
	  else contrib=PTPN*contrib0;
	  
	  (*out)=std::make_pair(neighbors[n],contrib);++out;	  
	}    

      return (1<<NDIM);
    }
    
    // Correct the gradient contribution to deltaWeight if necessary, when an edge is 
    // crossing a periodic boundary
    template <class T, class WF>
    Float getPeriodizedDeltaWeight(const T  * coord, const WF &wf, 
				   const IncidentEdge &edge)
    {
      Float result=edge.deltaWeight; //return result;
      
      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
	{
	  // We need to be carefull here, as the gradient contribution is computed
	  // relative to the vertex with index 0 reference, so this vertex has to
	  // be on the same side as coord even when the simplex crosses a periodic
	  // boundary ! If this is not the case, then we have to correct it !
	  const Coord* refC=edge.simplex->getVertex(0)->getCoordsConstPtr();
	  const Coord* otherRefC=refC;
	  
	  // otherSimplex may not exists on a boundary !
	  if (edge.otherSimplex!=NULL)
	    otherRefC=edge.otherSimplex->getVertex(0)->getCoordsConstPtr();

	  T ref[2];
	  T otherRef[2];	  
	  std::copy_n(refC,2,ref);
	  std::copy_n(otherRefC,2,otherRef);	 

	  T u0 = meshGeometry->checkCoordConsistency(ref[0],coord[0],0);
	  T u1 = meshGeometry->checkCoordConsistency(ref[1],coord[1],1);
	  T v0 = meshGeometry->checkCoordConsistency(otherRef[0],coord[0],0);
	  T v1 = meshGeometry->checkCoordConsistency(otherRef[1],coord[1],1);
	      
	  if ((u0!=ref[0])||(u1!=ref[1])||(v0!=otherRef[0])||(v1!=otherRef[1]))
	    {
	      // deltaWeight needs to be updated as this is where we store the 
	      // contribution of the reference point
	      /*
	      Float weight=wf.get(edge.simplex);	      	 
	      Float grad[2];
	      wf.getGradient(edge.simplex,grad);

	      // We compensate for the error we will do on the coordinate by first
	      // substracting it from the 0th order term.
	      weight -= 
		grad[0]*(u0-ref[0]) + 
		grad[1]*(u1-ref[1]);
	      */
	      
	      T coords[2]={u0,u1};
	      Float weight=wf.template compute<T,Float>(edge.simplex,coords);

	      Float otherWeight=0;
	      if (edge.otherSimplex!=NULL)
		{
		  // not a boundary edge
		  /*
		  otherWeight=wf.get(edge.otherSimplex);
		  Float otherGrad[2];	      
		  wf.getGradient(edge.otherSimplex,otherGrad);

		  otherWeight -= 
		    otherGrad[0]*(v0-otherRef[0]) + 
		    otherGrad[1]*(v1-otherRef[1]);
		  */
		  T coords[2]={v0,v1};
		  otherWeight=wf.template compute<T,Float>(edge.otherSimplex,coords);
		}
	  
	      if (edge.flags & IncidentEdge::flag_boundary) 
		otherWeight = 0; // boundary edge
	      else if (edge.flags & IncidentEdge::flag_fold) 
		otherWeight = -otherWeight; // fold edge
	      
	      result = weight-otherWeight;	      
	    }
	}
      return result;
    }

    // Correct the gradient contribution to 0th order weight if necessary, when a simplex
    // is crossing a periodic boundary
    template <class T, class WF>
    Float getPeriodizedWeight(const T  * coord, const WF &wf, 
			      Simplex *simplex, bool coordsAreConsistent=false)
    {
      Float result=wf.template get<Float>(simplex); //return result;

      if ((AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)&&(WF::ORDER>0))//&&(!coordsAreConsistent))
	{
	  // We need to be carefull here, as the gradient contribution is computed
	  // relative to the vertex with index 0 reference, so this vertex has to
	  // be on the same side as coord even when the simplex crosses a periodic
	  // boundary ! If this is not the case, then we have to correct it !
	  //const Coord* ref = simplex->getVertex(0)->getCoordsConstPtr();
	  T ref[NDIM];
	  for (int i=0;i<NDIM;++i)
	    ref[i]=simplex->getVertex(0)->getCoord(i);

	  T u0 = meshGeometry->checkCoordConsistency(ref[0],coord[0],0);
	  T u1 = meshGeometry->checkCoordConsistency(ref[1],coord[1],1);
		      
	  if ((u0!=ref[0])||(u1!=ref[1]))
	    {
	      T coords[2]={u0,u1};
	      result=wf.template compute<T,Float>(simplex,coords);
	      /*    	      
	      Float grad[2];
	      wf.getGradient(simplex,grad);
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
	      result -= 
		grad[0]*(u0-ref[0]) + 
		grad[1]*(u1-ref[1]);
	      */	      
	    }
	}
      return result;
    }


  };
  /*
  void getLeavesFromBBox(double simplexBBox[2][NDIM])
  {
    Voxel *origin=amr->getVoxelAt(simplexBBox[0]);
  }
  */
}

#include "../../internal/namespace.footer"
#endif
