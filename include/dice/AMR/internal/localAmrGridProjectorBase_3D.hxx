#ifndef __LOCAL_AMR_GRID_PROJECTOR_BASE_3D__HXX__
#define __LOCAL_AMR_GRID_PROJECTOR_BASE_3D__HXX__

#include <string>    
#include <iostream>   
#include <sstream>

#include "localAmrGridProjectorBasePrototype.hxx"

#include "../../internal/namespace.header"

//#define DEBUGDUMP 1
/*
#define DEBUG_CHECK 7971373578124785664
#define DEBUG_CHECK_RANK 0
*/
#ifdef DEBUGDUMP
#define PRINTVERTEX(v) {const Coord *crd=v->getCoordsConstPtr();fprintf(fl[0],"%e %e %e\n",crd[0],crd[1],crd[2]);}
#define PRINTCOORD(c) {fprintf(fl[0],"%e %e %e\n",c[0],c[1],c[2]);}
#endif

namespace internal {

  // unnamed namespace 
  namespace {
    // Generic version 
    template <int NDIM, class M, typename F, typename HF>
    struct FaceVectors;
    
    // Version where HF and F are the same
    template <class M, typename F>
    struct FaceVectors<3,M,F,F>
    {
      static const int NDIM=3;

      typedef F Float;
      typedef F HFloat;      

      typedef typename M::Simplex Simplex;
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord Coord;

      template <class Edge, class G>
      static void computeT(Edge &edge, G &meshGeometry, Float tolerance=1.E-3)
      {
	meshGeometry->template
	  getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
				      edge.otherVertex->getCoordsConstPtr(),
				      edge.vecT);
		  
	Float nrm = meshGeometry->template normalize_noCheck<Float,NDIM>(edge.vecT);	 
	if (nrm==0)
	  {
	    PRINT_SRC_INFO(LOG_WARNING);
	    glb::console->print<LOG_WARNING>("I encountered a highly degenerate tetrahedron (two vertices at the exact same location). Cross your fingers ;)\n");
	    edge.vecT[0]=1;edge.vecT[1]=0;edge.vecT[2]=0;
	  }
      }
      
      // Compute normal and binormal vectors for a given face adjacent to an edge
      template <class Face, class Edge, class G>
      static void computeNB(Face &face, Edge &edge, 
			    Vertex *faceVertex, Vertex *oppVertex, 
			    G &meshGeometry, Float tolerance=1.E-3)
      {
	Float dot;
	meshGeometry->template getVector<Coord,Float,NDIM>
	  (edge.vertex->getCoordsConstPtr(),
	   faceVertex->getCoordsConstPtr(),
	   face.vecN);
      
	dot = meshGeometry->template dot_noCheck<Float,Float,NDIM>(edge.vecT,face.vecN);
	for (int i=0;i<NDIM;++i) face.vecN[i]-=dot*edge.vecT[i];
	dot = meshGeometry->template normalize_noCheck<Float,NDIM>(face.vecN);

	// And cross product for vecB (NOT gram-schmidt as it would break when a simplex
	// is almost or exactly flat)      
	face.vecB[0]=edge.vecT[1]*face.vecN[2]-edge.vecT[2]*face.vecN[1];
	face.vecB[1]=edge.vecT[2]*face.vecN[0]-edge.vecT[0]*face.vecN[2];
	face.vecB[2]=edge.vecT[0]*face.vecN[1]-edge.vecT[1]*face.vecN[0];          

	// now orient face.vecB toward the remaining vertex in the simplex: the one that 
	// does not belong to the face (i.e. oppVertex)
	Float vecB[NDIM];
	meshGeometry->template getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
							   oppVertex->getCoordsConstPtr(),
							   vecB);
	dot=meshGeometry->template dot_noCheck<Float,Float,NDIM>(face.vecB,vecB);
	if (dot<0)
	  {
	    for (int i=0;i<NDIM;++i) 
	      face.vecB[i]=-face.vecB[i];
	  }      
      }
    };

    // F and HF are different floating point types with low and high precision respectively.
    // HF is used only if F is not precise enough.
    template <class M, typename F, typename HF>
    struct FaceVectors<3,M,F,HF>
    {
      static const int NDIM=3;

      typedef F Float;
      typedef HF HFloat;      

      typedef typename M::Simplex Simplex;
      typedef typename M::Vertex  Vertex;
      typedef typename M::Coord Coord;

      template <class Edge, class G>
      static void computeT(Edge &edge, G *meshGeometry, Float tolerance=1.E-3)
      {
	FaceVectors<NDIM,M,Float,Float>::
	  computeT(edge,meshGeometry);

	bool recomp=false;
	for (int i=0;i<NDIM;++i)
	  if (fabs(edge.vecT[i]) < tolerance)
	    recomp=true;

	if (recomp) 
	  {
	    HFloat vecHT[NDIM];
	    meshGeometry->template
	      getVector<Coord,HFloat,NDIM>(edge.vertex->getCoordsConstPtr(),
					   edge.otherVertex->getCoordsConstPtr(),
					   vecHT);
	    HFloat nrm=meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHT);
	    if (nrm==0)
	      {
		PRINT_SRC_INFO(LOG_WARNING);
		glb::console->print<LOG_WARNING>("I encountered a highly degenerate tetrahedron (two vertices at the exact same location). Cross your fingers ;)\n");
		edge.vecT[0]=1;edge.vecT[1]=0;edge.vecT[2]=0;
	      }

	    for (int i=0;i<NDIM;++i)
	      edge.vecT[i]=hlp::numericStaticCast<Float>(vecHT[i]);
	  }	
      }
      
      template <class Face, class Edge, class G>
      static void computeNB(Face &face, Edge &edge, 
			    Vertex *faceVertex, Vertex *oppVertex, 
			    G *meshGeometry, Float tolerance=1.E-3)
      {
	FaceVectors<NDIM,M,Float,Float>::
	  computeNB(face,edge,faceVertex,oppVertex,meshGeometry);

	bool recomp=false;
	for (int i=0;i<NDIM;++i)
	  if ((fabs(face.vecB[i]) < tolerance)||
	      (fabs(face.vecN[i]) < tolerance))
	    recomp=true;

	if (recomp) 
	  {
	    HFloat dotH;
	    HFloat vecHN[NDIM];
	    HFloat vecHT[NDIM];
	    
	    // We always recompute vecHT as we need it at HFloat precision.
	    meshGeometry->template
	      getVector<Coord,HFloat,NDIM>(edge.vertex->getCoordsConstPtr(),
					   edge.otherVertex->getCoordsConstPtr(),
					   vecHT);
		  
	    meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHT);

	    /*
	    for (int i=0;i<NDIM;++i)
	      vecHT[i]=hlp::numericStaticCast<HFloat>(edge.vecT[i]);
	    */

	    // Gram-Schmidt for vecN
	    // FIXME: this may be a problem if a vertex falls on a segment ???
	    meshGeometry->template getVector<Coord,HFloat,NDIM>
	      (edge.vertex->getCoordsConstPtr(),
	       faceVertex->getCoordsConstPtr(),
	       vecHN);  

	    dotH = meshGeometry->template dot_noCheck<HFloat,HFloat,NDIM>(vecHT,vecHN);
	    for (int i=0;i<NDIM;++i) vecHN[i]-=dotH*vecHT[i];
	    dotH = meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHN);

	    // And cross product for vecB (NOT gram-schmidt as it would break when a simplex
	    // is almost or exactly flat)
	    HFloat vecHB[NDIM];
	    vecHB[0]=vecHT[1]*vecHN[2]-vecHT[2]*vecHN[1];
	    vecHB[1]=vecHT[2]*vecHN[0]-vecHT[0]*vecHN[2];
	    vecHB[2]=vecHT[0]*vecHN[1]-vecHT[1]*vecHN[0];
      
	    // now orient face.vecB toward the remaining vertex in the simplex: 
	    // the one that does not belong to the face (i.e. oppVertex)
	    HFloat vecHB_dir[NDIM];
	    meshGeometry->template getVector<Coord,HFloat,NDIM>
	      (edge.vertex->getCoordsConstPtr(),
	       oppVertex->getCoordsConstPtr(),
	       vecHB_dir);
	    
	    dotH=meshGeometry->template dot_noCheck<HFloat,HFloat,NDIM>(vecHB,vecHB_dir);

	    if (dotH<0)
	      {
		for (int i=0;i<NDIM;++i)
		  {
		    face.vecB[i]=-hlp::numericStaticCast<Float>(vecHB[i]);
		    face.vecN[i]=hlp::numericStaticCast<Float>(vecHN[i]);
		  }
	      }
	    else
	      {
		for (int i=0;i<NDIM;++i)
		  {
		    face.vecB[i]=hlp::numericStaticCast<Float>(vecHB[i]);
		    face.vecN[i]=hlp::numericStaticCast<Float>(vecHN[i]);
		  }
	      }
	  }
	
      }
    };
  }


  template <class AMR, class MESH, typename F, typename HF, bool checkAccuracy, int dummy> 
  class LocalAmrGridProjectorBaseT<3,AMR,MESH,F,HF,checkAccuracy,dummy>
  {
  public:
    static const int NDIM = 3;

    typedef typename MESH::vertexPtr_iterator vertexPtr_iterator;
    typedef typename MESH::simplexPtr_iterator simplexPtr_iterator;
    typedef typename MESH::Simplex Simplex;
    typedef typename MESH::Vertex  Vertex;
    typedef typename MESH::Segment Segment;
    typedef typename MESH::SegmentHandle SegmentHandle;
    typedef typename MESH::FacetHandle FacetHandle;
    typedef typename MESH::segment_circulator segment_circulator;
    typedef typename MESH::Coord Coord;
    typedef typename MESH::GeometricProperties MeshGeometricProperties;

    typedef typename AMR::ICoord ICoord;
    typedef typename AMR::Voxel Voxel;
    typedef typename AMR::GeometricProperties AmrGeometricProperties;

    typedef F Float;
    typedef HF HFloat;

    LocalAmrGridProjectorBaseT(AMR *amr_, MESH *mesh_, 
			       int nThreads_)://,MpiCommunication *mpiCom_):
      nThreads(nThreads_),
      //mpiCom(mpiCom_),
      amr(amr_),
      mesh(mesh_)
    {
#ifdef DEBUG_CHECK
      mpiCom=glb::mpiComWorld;
#endif
      meshGeometry = mesh->getGeometry();
      amrGeometry = amr->getGeometry();    
#ifdef DEBUGDUMP
      nFiles=0;      
      openFile("vectors");
#endif
    }

    virtual ~LocalAmrGridProjectorBaseT()
    {
#ifdef DEBUGDUMP
      closeFile();
#endif
    }   

  protected:
    int nThreads;
#ifdef DEBUG_CHECK
    MpiCommunication *mpiCom;
#endif
    AMR *amr;
    MESH *mesh;
    MeshGeometricProperties *meshGeometry;
    AmrGeometricProperties *amrGeometry;    
    
  public:
    struct IncidentFace {
      static const int flag_boundary = (1<<0);
      static const int flag_invalidNext = (1<<1);
      static const int flag_invalidPrev = (1<<2);
      //static const int flag_boundary2 = (1<<1);
      //static const int flag_boundary = ((1<<0)|(1<<1));
      static const int flag_fold = (1<<3);

      Simplex *prevSimplex; // the simplex in -vecB direction 
      Simplex *nextSimplex; // the simplex in vecB direction
      Float deltaWeight;
      Float deltaGrad[NDIM];   
      //Float prevWeight;
      //Float nextWeight;
      //Float foldingFactor;
      Float vecN[NDIM];     // Normal vector
      Float vecB[NDIM];     // Binormal vector
      int flags;
    };

    struct IncidentEdge {
      Vertex *vertex;
      Vertex *otherVertex;
      //Float otherContrib; // contribution to otherVertex
      Simplex *owner;      // useful ?

      Float vecT[NDIM];    // Tangent vector      
      std::vector<IncidentFace> face; // all the simplex facets incident to the edge

      template <class L>
      void print() const
      {
	std::stringstream ss;
	// for (int i=0;i<face.size();++i) 
	//   ss<<"["<<face[i].prevWeight<<","<<face[i].nextWeight<<"]";	
	// for (int i=0;i<face.size();++i) 
	//   ss<<"["<<face[i].deltaWeight<<"]";

	for (int i=0;i<face.size();++i) 
	  ss<<"["<<face[i].deltaWeight<<",("<<
	    face[i].deltaGrad[0]<<","<<
	    face[i].deltaGrad[1]<<","<<
	    face[i].deltaGrad[2]<<"),("<<	    
	    face[i].vecN[0]<<","<<
	    face[i].vecN[1]<<","<<
	    face[i].vecN[2]<<"),"<<"),("<<	    
	    face[i].vecB[0]<<","<<
	    face[i].vecB[1]<<","<<
	    face[i].vecB[2]<<"),"<<
	    face[i].flags<<"]\n";
	    
	

	glb::console->print<L>("EDGE:%ld facets:%s\n",face.size(),ss.str().c_str());
      }
    };
    
    template <class WF>
    void getVertexContribAndEdges(Simplex * const * simplex, int nSimplices,
				   Vertex *vertex, 
				   const WF &wf,
				   std::vector<IncidentEdge> &ownedEdges,
				   bool debug=false)
    {    
      ownedEdges.clear();
      // Recover the edges we own
      for (int i=0;i<nSimplices;++i)
	getEdges(simplex[i],vertex,wf,ownedEdges);

      // and compute the contribution to vertex, as well as to the other vertex in each
      // owned segment
      /*
      Float vertexContrib=0;  
      for (unsigned long i=0;i<ownedEdges.size();++i)
	{
	  Float maxLocalContrib;
	  Float contrib=computeContrib(vertex->getCoordsConstPtr(),
				       wf,ownedEdges[i],maxLocalContrib);	  
	  // contrib for the opposite vertex
	  //ownedEdges[i].otherContrib = -contrib; 
	  vertexContrib += contrib;
	  if (maxContrib<maxLocalContrib) 
	    maxContrib=maxLocalContrib;
	}

       if (debug)
	{
	  Float c=0;
	  //std::cout.setf( std::ios::fixed );
	  std::cout.precision(15);
	  std::cout << std::scientific;	  
	  for (int i=0;i<ownedEdges.size();++i)
	    {
	      ownedEdges[i].template print<LOG_STD_ALL>();
	      Float maxLocalContrib;
	      Float contrib=
		computeContrib(vertex->getCoordsConstPtr(),
			       wf,ownedEdges[i],maxContrib,true);
	      
	      
	      std::cout << "c = " << c  << " + " << contrib << " = " << contrib+c <<"\n";
	      std::cout << "Contrib = " << contrib 
			<< "/ Max contrib = " << maxContrib 
			<< ((contrib<maxLocalContrib*1.E-4)?"OKOKOK":"WRONG!")<<"\n";
	      c+=contrib;	      
	    }
	  
	}             
      
      return vertexContrib;   
      */
      //return 0;
    }
    
    template <class WF>
    void getEdges(Simplex *simplex, Vertex *vertex, 
		  const WF &wf, std::vector<IncidentEdge> &ownedEdges)
    {
      if (invalidSimplex(simplex)) return;

      for (int i=0;i<Simplex::NVERT;++i)
	{
	  // First check which edges of the simplex departing from vertex we own (if any)

	  // Condition1: the current vertex's adress is higher than that of the vertex
	  // at the other extremity
	  if (simplex->getVertex(i)<vertex)
	    {	      
	      segment_circulator ci=
		Segment(vertex,simplex->getVertex(i),simplex).getCirculator();	    
	      
	      char tag=wf.getTag(simplex);	      
	      int nTagged=0;
	      int nFaces=0;

	      if (!tag)
		{
		  do {++ci;++nFaces;} 
		  while ((*ci < simplex)||invalidSimplex(*ci));
		}
	      else
		{
		  do {
		    ++ci;
		    ++nFaces;
		    if ((!(*ci)->isLocal())||(wf.getTag(*ci))) ++nTagged;
		  } 
		  while ((*ci < simplex)||invalidSimplex(*ci));
		}
	      //while ((*ci < simplex)||(ci->isGhost()));
	      

	      /*
		// This works when we don't remove the tagged simplices
	      int nFaces=0;
	      int nTagged=0;
	      do {++ci;nFaces++;} 
	      while ((*ci < simplex)||invalidSimplex(*ci));
	      */

	      // Condition2: the current simplex has an address higher than any other 
	      // simplex around the current segment (-> we made a full turn around the
	      // edge with the segment circulator and came back to simplex)
	      if ((*ci == simplex)&&(nFaces!=nTagged))
		{
		  ownedEdges.resize(ownedEdges.size()+1);
		  IncidentEdge &curEdge = ownedEdges.back();
		  curEdge.vertex = vertex;
		  curEdge.otherVertex = simplex->getVertex(i);
		  curEdge.owner = simplex;
		  curEdge.face.resize(nFaces);
		  nFaces = 0;
		  //for (long i=0;i<curEdge.face.size();++i)
		  //curEdge.face[i].flags=0;

		  // printf("ANDNET\n3\n%d\n",nFaces*4+nFaces*4);

		  FaceVectors<NDIM,MESH,Float,HFloat>::computeT(curEdge,meshGeometry);
		  /* // This was replaced by computeT
		  meshGeometry->template
		    getVector<Coord,Float,NDIM>(vertex->getCoordsConstPtr(),
						curEdge.otherVertex->getCoordsConstPtr(),
						curEdge.vecT);
		  
		  Float nrm = meshGeometry->template normalize_noCheck<Float,NDIM>
		    (curEdge.vecT);
		  
		  if (nrm==0)
		    {
		      PRINT_SRC_INFO(LOG_WARNING);
		      glb::console->print<LOG_WARNING>("I encountered a highly degenerate tetrahedron (two vertices at the exact same location). Cross your fingers ;)\n");
		      curEdge.vecT[0]=1;curEdge.vecT[1]=0;curEdge.vecT[2]=0;
		    }
		  */

		  // Now compute the normal and binormal to each facet around the edge
		  Simplex *oldSimplex;
		  Simplex *newSimplex=simplex;
		  //ci=ci_end;		  
		  do
		    {
		      oldSimplex=newSimplex;
		      int oldFacetIndex=ci.getNextFacetIndex();
		      ++ci;newSimplex=*ci;

		      if (ci.crossedBoundary())
			{			  
			  // We crossed a boundary, so we have 2 distinct contributions
			  // to compute (1 per boundary face) 
			  // => we must add an extra contribution
			  curEdge.face.resize(curEdge.face.size()+1); 
			  //curEdge.face.back().flags=0;

			  Vertex *faceVertex=NULL;
			  Vertex *oppVertex = oldSimplex->getVertex(oldFacetIndex);
			  //int facetId=oldFacetIndex;
			  for (int j=0;j<Simplex::NVERT;++j)
			    {
			      Vertex *curVertex=oldSimplex->getVertex(j);
			      
			      if ((curVertex!=oppVertex)&&
				  (curVertex!=vertex)&&
				  (curVertex!=curEdge.otherVertex))
				faceVertex = curVertex;
			    }
			  //curEdge.face[nFaces].flags |= IncidentFace::flag_boundary1;
			  getFace(curEdge,oldSimplex,NULL,faceVertex,
				  oppVertex,wf,curEdge.face[nFaces++]);
			  

			  faceVertex = newSimplex->getVertex(ci.getNextFacetIndex());

			  for (int j=0;j<Simplex::NVERT;++j)
			    {
			      Vertex *curVertex=newSimplex->getVertex(j);
			      if ((curVertex!=faceVertex)&&
				  (curVertex!=vertex)&&
				  (curVertex!=curEdge.otherVertex))
				{
				  //facetId = j;
				  oppVertex=curVertex;
				}
			    }
			  //curEdge.face[nFaces].flags |= IncidentFace::flag_boundary2;
			  getFace(curEdge,newSimplex,NULL,faceVertex,
				  oppVertex,wf,curEdge.face[nFaces++]);
			}
		      else
			{
			  
			  if (newSimplex->isLocal()&&
			      oldSimplex->isLocal()&&
			      wf.getTag(newSimplex)&&
			      wf.getTag(oldSimplex))
			    {
			      // The two simplices are tagged, so we can remove that face
			      curEdge.face.resize(curEdge.face.size()-1); 
			    }
			  else
			    {
			      //int facetId;
			      Vertex *oppVertex=NULL;
			      Vertex *faceVertex=NULL;
			      for (int j=0;j<Simplex::NVERT;++j)
				{
				  Vertex *curVertex=newSimplex->getVertex(j);
				  if (newSimplex->getNeighbor(j)==oldSimplex) 
				    {
				      oppVertex = curVertex;
				      //facetId = j;
				    }
				  else if ((curVertex!=vertex)&&
					   (curVertex!=curEdge.otherVertex))
				    faceVertex=curVertex;
				}

			      getFace(curEdge,newSimplex,oldSimplex,faceVertex,
				      oppVertex,wf,curEdge.face[nFaces]);
			      nFaces++;
			    }
			    /*
			  for (int j=0;j<Simplex::NVERT;++j)
			    PRINTVERTEX(newSimplex->getVertex(j));
			  Coord vc[Vertex::NDIM_W];
			  Float fac=1.0/128.0;
			  curEdge.vertex->getCoords(vc);
			  PRINTCOORD(vc);
			  for (int j=0;j<NDIM;++j) vc[j]+=fac*curEdge.vecT[j];
			  PRINTCOORD(vc);
			  for (int j=0;j<NDIM;++j) vc[j]+=fac*curEdge.face[nFaces-1].vecN[j];
			  PRINTCOORD(vc);
			  for (int j=0;j<NDIM;++j) vc[j]+=fac*curEdge.face[nFaces-1].vecB[j];
			  PRINTCOORD(vc);
			  */
			}

		      //} while (ci != ci_end);
		    } while (*ci != simplex);
		
		  /*		
		  printf("3 %d\n",nFaces);
		  for (int j=0;j<nFaces;++j) 
		    printf("%d %d %d %d\n",8*j,8*j+1,8*j+2,8*j+3);
		  printf("1 %d\n",nFaces*3);
		  for (int j=0;j<nFaces;++j) 
		    printf("%d %d\n%d %d\n%d %d\n",8*j+4,8*j+5,8*j+5,8*j+6,8*j+6,8*j+7);
		  */
		} // We own the segment
	      
	    }
	}
      
    }

    template <typename R, class WF, class OutputIterator>
    int addVoxelFacetSimplexEdgeContrib(const R &raytracer, const WF &wf,
					const IncidentEdge &edge, 
					OutputIterator out)
    {  
      typedef typename R::Coord RCoord;

      const RCoord *exitPoint  = raytracer.getExitPointConstPtr();
      Voxel  *curVoxel   = raytracer.getCurVoxel();
      Voxel  *nextVoxel  = raytracer.getNextVoxel();
      int     normalDim  = raytracer.getExitFaceNormalDim();
      Float  normalSign = raytracer.getExitFaceNormalSign();

      int dbg=0;
#ifdef DEBUG_CHECK
      const ICoord ic=DEBUG_CHECK;
      int checkRank=DEBUG_CHECK_RANK;
      bool checkNext=false;
      bool checkCur=false;
      if ((curVoxel->getIndex()==ic)&&(mpiCom->rank()==checkRank))
	checkCur=true;
      if ((nextVoxel->getIndex()==ic)&&(mpiCom->rank()==checkRank))
	checkNext=true;

      if (checkCur||checkNext) {dbg=1;edge.template print<LOG_STD_ALL>();}
#endif
      
      /*
      if ((edge.vertex->getLocalIndex()==371438)||
	  (edge.vertex->getLocalIndex()==177696))
	dbg++;
      if ((edge.otherVertex->getLocalIndex()==371438)||
	  (edge.otherVertex->getLocalIndex()==177696))
	dbg++;
      if (dbg<2) dbg=0;
      if (dbg) {std::cout.precision(20);std::cout << std::scientific;}
      */
      Float maxContrib=0;//used to store the maximal contrib, to control accuracy
      if (AMR::BOUNDARY_TYPE != BoundaryType::PERIODIC)
	{
	  Float contrib = 
	    computeVoxelFacetSimplexEdgeContrib
	    (exitPoint,normalDim,normalSign,wf,edge,
	     maxContrib,dbg);

	  /*DELETME*/
	  /*
	  if (fabs(fabs(contrib)-6.768231E-4)<1.E-6)
	    {
	      Float c2=computeVoxelFacetSimplexEdgeContrib
	      (exitPoint,normalDim,normalSign,wf,edge,true);
	      printf("(%e==%e) -> %d\n",contrib,c2,contrib==c2);
	    }
	    */
	  if (checkAccuracy)
	    {
	      contrib-=maxContrib;	      
	      (*out) = std::make_pair(curVoxel,-maxContrib);++out;
	      (*out) = std::make_pair(nextVoxel,maxContrib);++out;
	    }
	  (*out) = std::make_pair(curVoxel,-contrib);++out;
	  (*out) = std::make_pair(nextVoxel,contrib);++out;
	  
#ifdef DEBUG_CHECK
	  if (checkCur) 
	    std::cout << "ADDING ray curVoxel contrib : " << -contrib << std::endl;
	  if (checkNext) 
	    std::cout << "ADDING ray newVoxel contrib : " << contrib << std::endl;	    
#endif
	}
      else
	{
	  int boundary=0;

	  for (int i=0;i<NDIM;++i) 
	    {
	      if (amr->getBBoxMax(i)-amr->getEpsilon(i)<exitPoint[i]) 
		boundary|=1<<i;
	      else if (amr->getBBoxMin(i)+amr->getEpsilon(i)>exitPoint[i]) 
		boundary|=1<<i;
	    }

	  if (boundary)
	    {
	      RCoord refPoint[NDIM];
	      RCoord exitPoint_Cur[NDIM];
	      RCoord exitPoint_Next[NDIM];

	      amr->index2CenterCoords(curVoxel->getIndex(),refPoint);
	      for (int i=0;i<NDIM;++i)
		exitPoint_Cur[i]= amrGeometry->template 
		  checkCoordConsistency<RCoord>(exitPoint[i],refPoint[i],i);
		  // refPoint[i]+ amrGeometry->template 
		  // correctCoordsDiff<RCoord>(exitPoint[i]-refPoint[i],i);

	      amr->index2CenterCoords(nextVoxel->getIndex(),refPoint);
	      for (int i=0;i<NDIM;++i)
		exitPoint_Next[i]=amrGeometry->template 
		  checkCoordConsistency<RCoord>(exitPoint[i],refPoint[i],i);
		  // refPoint[i]+ amrGeometry->template 
		  // correctCoordsDiff<RCoord>(exitPoint[i]-refPoint[i],i);

	      /*
	      if (findVertex(0.328834,0.238395,-1,exitPoint))
		{
		  computeVoxelFacetSimplexEdgeContrib
		    (exitPoint_Cur,normalDim,normalSign,wf,edge,true);
		}
	      */
	      //if (dbg) std::cout << "On Boundary!"<<std::endl;
	      Float contribCur = 
		computeVoxelFacetSimplexEdgeContrib
		(exitPoint_Cur,normalDim,normalSign,wf,edge,maxContrib,dbg);

	      if (checkAccuracy)
		{
		  contribCur-=maxContrib;
		  (*out)=std::make_pair(curVoxel,-maxContrib);++out;
		}
	      (*out)=std::make_pair(curVoxel,-contribCur);++out;
	      

	      Float contribNext = 
		computeVoxelFacetSimplexEdgeContrib
		(exitPoint_Next,normalDim,normalSign,wf,edge,maxContrib,dbg);
	      
	      if (checkAccuracy)
		{
		  contribNext-=maxContrib;
		  (*out)=std::make_pair(nextVoxel,maxContrib);++out;
		}
	      (*out)=std::make_pair(nextVoxel,contribNext);++out;

#ifdef DEBUG_CHECK
	      if (checkCur) 
		std::cout << "ADDING ray curVoxel contrib (pb): " << -contribCur 
			  << std::endl;
	      if (checkNext) 
		std::cout << "ADDING ray newVoxel contrib (pb): " << contribNext 
			  << std::endl;	  
#endif

	      // if (boundary&(1<<0))
	      // printf("1:(%g,%g,%g) -> %g\n",exitPoint_Cur[0],exitPoint_Cur[1],exitPoint_Cur[2],-contribCur);
	      // printf("2:(%g,%g,%g) -> %g\n",exitPoint_Next[0],exitPoint_Next[1],exitPoint_Next[2],contribNext);
	    }
	  else
	    {
	      Float contrib = computeVoxelFacetSimplexEdgeContrib
		(exitPoint,normalDim,normalSign,wf,edge,maxContrib,dbg);   
	      
	      if (checkAccuracy)
		{
		  contrib-=maxContrib;
		  (*out) = std::make_pair(curVoxel,-maxContrib);++out;
		  (*out) = std::make_pair(nextVoxel,maxContrib);++out;
		}
	      (*out) = std::make_pair(curVoxel,-contrib);++out;
	      (*out) = std::make_pair(nextVoxel,contrib);++out;

	      

#ifdef DEBUG_CHECK
	      if (checkCur) 
		std::cout << "ADDING ray curVoxel contrib (p): " << -contrib << std::endl;
	      if (checkNext) 
		std::cout << "ADDING ray nextVoxel contrib (p): " << contrib << std::endl;
#endif
	    }
	}
      
      return 4;
    }
    
    // projectedFacetNormal is the normalized normal to the facet. Its direction
    // is not relevant (it can be pointing inward or outward).
    // vEdgeSign must be either +1 or -1
    // NB: the voxel edge orientation (in contribSign[][dim]) is modified
    // by this function so that it becomes aligned with the facetNormal
    // NB2: the order of segNeighbors must respect the convention of LocalAMRgrid !
    template <class WF, class OutputIterator>
    int addVoxelEdgeSimplexFacetContrib(Simplex *simplex, 
					int simplexNeighborIndex, 
					Float projectedFacetNormal[NDIM], 
					int dim,
					Float vEdgeSign[NDIM], 
					Float intersection[NDIM],
					Voxel *segNeighbor[Voxel::NSEGNEI],	
					const WF &wf, OutputIterator out,
					bool coordsAreConsistent=false,
					bool dbg=false) const
    {   
      //int nContribs=5;
      Float maxContrib=0;
      Float contrib=computeVoxelEdgeSimplexFacetContrib
	(simplex,simplexNeighborIndex,projectedFacetNormal,dim,
	 vEdgeSign,intersection,wf,maxContrib,coordsAreConsistent,dbg);
      
      if (checkAccuracy)
	{
	  (*out)=std::make_pair(segNeighbor[0],contrib-maxContrib);++out;
	  (*out)=std::make_pair(segNeighbor[0],maxContrib);++out;
	}
      else
	{
	  (*out)=std::make_pair(segNeighbor[0],contrib);++out;
	}
    
#ifdef DEBUG_CHECK
      const ICoord ic=DEBUG_CHECK;
      int checkRank=DEBUG_CHECK_RANK;
      int check=0;
      for (int i=0;i<4;++i) 
	if ((segNeighbor[i]->getIndex()==ic)&&(mpiCom->rank()==checkRank)) 
	  {
	    check=i+1;
	    Simplex *nei = simplex->getNeighbor(simplexNeighborIndex);
	    Float wOut=validSimplex(nei)?nei->cache.d:0;
	    //Float wIn=simplex->cache.d;
	    Float wIn=validSimplex(simplex)?simplex->cache.d:0; 

	    std::cout << "Voxel edge contrib : wOut="<< wOut
		      << ", wIn="<< wIn
		      << " ("<<intersection[0]<<","<<intersection[1]<<","<<intersection[2]
		      << ")" << std::endl;

	    std::cout << "Proj normal = ["<<projectedFacetNormal[0]<<","
		      <<projectedFacetNormal[1]<<","
		      <<projectedFacetNormal[2]<<"]"<<std::endl;
	    // printf("Voxel edge contrib : wOut=%Lg, wIn=%Lg (%lg,%lg,%lg)\n",wOut,wIn,
	    // 	   intersection[0],intersection[1],intersection[2]);
	    // printf("proj normal = [%lg %lg %lg]\n",
	    // 	   projectedFacetNormal[0],
	    // 	   projectedFacetNormal[1],
	    // 	   projectedFacetNormal[2]);
	    /*
	    contrib=computeVoxelEdgeSimplexFacetContrib
	      (simplex,simplexNeighborIndex,projectedFacetNormal,
	       dim,vEdgeSign,intersection,true);
	    */
	  } 
#endif     

      // if boundary conditions are periodic, the intersection coordinates have to change
      // depending on the segNeighbor when the segment is along a boundary ...
      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
	{
	  int boundaryIndex[NDIM];	 
	  int alongBoundary=0;
	  for (int d=0;d<NDIM;++d)
	    {	      
	      if (d!=dim)
		{
		  Float c=intersection[d];
		  if ((c>=amr->getBBoxMax(d)-amr->getEpsilon(d))||
		      (c<=amr->getBBoxMin(d)+amr->getEpsilon(d)))
		    {
		      boundaryIndex[alongBoundary++]=d;		     
		    }
		}
	    }

	  if (alongBoundary==0)
	    {
	      // Nothing special to check !
	      // (*out)=std::make_pair(segNeighbor[1],-contrib);++out;
	      // (*out)=std::make_pair(segNeighbor[2],-contrib);++out;	  
	      // (*out)=std::make_pair(segNeighbor[3],contrib);++out;
	      if (checkAccuracy)
		{
		  contrib-=maxContrib;
		  (*out)=std::make_pair(segNeighbor[1],-maxContrib);++out;
		  (*out)=std::make_pair(segNeighbor[2],-maxContrib);++out;	  
		  (*out)=std::make_pair(segNeighbor[3],maxContrib);++out;
		}
	      (*out)=std::make_pair(segNeighbor[1],-contrib);++out;
	      (*out)=std::make_pair(segNeighbor[2],-contrib);++out;	  
	      (*out)=std::make_pair(segNeighbor[3],contrib);++out;

	      
	    }
	  else 
	    {
	      Float maxPContrib[3];
	      Float pContrib[3];
	      //Float projectedFacetNormal2[NDIM];
	      Float intersection2[NDIM];

	      //std::copy_n(projectedFacetNormal,NDIM,projectedFacetNormal2);
	      
	      if (alongBoundary==1)
		{
		  // We are crossing one boundary, we can spare few computation ...
		  int bi=boundaryIndex[0]; // the dimension we are crossing
		  int otherDim=0; // the dimension we are not crossing

		  if (otherDim==bi) otherDim++;
		  if (otherDim==dim) 
		    {
		      otherDim++;
		      if (otherDim==bi) otherDim++;
		    }		      

		  std::copy(intersection,intersection+NDIM,intersection2);
		  intersection2[bi]+=amrGeometry->getBBoxSize(bi)*((vEdgeSign[bi]>0)?1:-1);
		  
		  Float contrib2=computeVoxelEdgeSimplexFacetContrib
		    (simplex,simplexNeighborIndex,projectedFacetNormal,dim,
		     vEdgeSign,intersection2,wf,maxPContrib[0]);
	      
		  // Because that is the order of the neighbors by convention.
		  if (bi<otherDim) 
		    {		     		      
		      maxPContrib[1]=maxContrib;
		      maxPContrib[2]=maxPContrib[0];
		      pContrib[0]=contrib2;
		      pContrib[1]=contrib;
		      pContrib[2]=contrib2;
		    }
		  else
		    {
		      maxPContrib[1]=maxPContrib[0];
		      maxPContrib[2]=maxPContrib[0];
		      maxPContrib[0]=maxContrib;
		      pContrib[0]=contrib;
		      pContrib[1]=contrib2;
		      pContrib[2]=contrib2;
		    }	
		}
	      else
		{
		  // We are crossing two boundaries
		  // and boundaryIndex[0]<boundaryIndex[1] by contruction
		  int b1=boundaryIndex[0];
		  int b2=boundaryIndex[1];
		  for (int i=1;i<(1<<2);++i)
		    {
		      std::copy(intersection,intersection+NDIM,intersection2);
		      if (i&(1<<0))
			{
			  intersection2[b1]+=
			    amrGeometry->getBBoxSize(b1) * ((vEdgeSign[b1]>0)?1:-1);
			  
			}
		      if (i&(1<<1))
			{
			  intersection2[b2]+=
			    amrGeometry->getBBoxSize(b2) * ((vEdgeSign[b2]>0)?1:-1);
			    
			}
		      pContrib[i-1]=computeVoxelEdgeSimplexFacetContrib
			(simplex,simplexNeighborIndex,projectedFacetNormal,dim,
			 vEdgeSign,intersection2,wf,maxPContrib[i-1]);		     
		    }
		}

	      if (checkAccuracy)
		{
		  pContrib[0]-=maxPContrib[0];
		  pContrib[1]-=maxPContrib[1];
		  pContrib[2]-=maxPContrib[2];
		  
		  (*out)=std::make_pair(segNeighbor[1],-maxPContrib[0]);++out;
		  (*out)=std::make_pair(segNeighbor[2],-maxPContrib[1]);++out;	  
		  (*out)=std::make_pair(segNeighbor[3],maxPContrib[2]);++out;
		}
	      (*out)=std::make_pair(segNeighbor[1],-pContrib[0]);++out;
	      (*out)=std::make_pair(segNeighbor[2],-pContrib[1]);++out;	  
	      (*out)=std::make_pair(segNeighbor[3],pContrib[2]);++out;
	      
	      //nContribs+=3;
	    }	   
	}
      else
	{
	  if (checkAccuracy)
	    {
	      contrib-=maxContrib;
	      
	      (*out)=std::make_pair(segNeighbor[1],-maxContrib);++out;
	      (*out)=std::make_pair(segNeighbor[2],-maxContrib);++out;	  
	      (*out)=std::make_pair(segNeighbor[3],maxContrib);++out;
	    }
	  (*out)=std::make_pair(segNeighbor[1],-contrib);++out;
	  (*out)=std::make_pair(segNeighbor[2],-contrib);++out;	  
	  (*out)=std::make_pair(segNeighbor[3],contrib);++out;	  

#ifdef DEBUG_CHECK
	  if (check)
	    {
	      if ((check==1)||(check==4))
		std::cout<<"ADDING V edge/S face : "<<contrib<<std::endl;
	      else
		std::cout<<"ADDING V edge/S face : "<<-contrib<<std::endl;	   
	    }
#endif
	}
      
      return 8;
    }

    template <int sign, class CT, class WF, class OutputIterator>
    int addContrib(const CT  * coord, const WF &wf, const IncidentEdge &edge,
		   Voxel *v, OutputIterator out, bool debug=false)
    {
      Float maxContrib;
      Float contrib=computeContrib(coord,wf,edge,maxContrib,debug);

      // We separate max contrib for accuracy checking
      if (checkAccuracy) contrib-=maxContrib;
      if (sign>0)
	{
	  (*out)=std::make_pair(v,contrib);++out;
	  if (checkAccuracy) {(*out)=std::make_pair(v,maxContrib);++out;}
	}
      else
	{
	  (*out)=std::make_pair(v,-contrib);++out;
	  if (checkAccuracy) {(*out)=std::make_pair(v,-maxContrib);++out;}
	}
      return 2;
    }

    // Compute the 2*nFaces contribs from a point ray intersection    
    template <class CT, class WF>
    Float computeContrib(const CT  * coord, const WF &wf, const IncidentEdge &edge,
			 Float &maxContrib,bool debug=false)
    {             
      maxContrib=0;
      Float contrib=0;
      Float PT= 
	coord[0] * edge.vecT[0]+
	coord[1] * edge.vecT[1]+
	coord[2] * edge.vecT[2];

      maxContrib=0;
    
      for (int f=0;f<edge.face.size();++f)
	{
	  const IncidentFace &curFace=edge.face[f];
	  
	  Float PN= 
	    coord[0] * curFace.vecN[0]+
	    coord[1] * curFace.vecN[1]+
	    coord[2] * curFace.vecN[2];

	  Float PB=
	    coord[0] * curFace.vecB[0]+
	    coord[1] * curFace.vecB[1]+
	    coord[2] * curFace.vecB[2];

	  Float deltaWeight = curFace.deltaWeight;	  
	  
	  if (WF::ORDER>0)
	    {
	      Float GT=
		edge.vecT[0]*curFace.deltaGrad[0]+
		edge.vecT[1]*curFace.deltaGrad[1]+
		edge.vecT[2]*curFace.deltaGrad[2];
	      Float GN=
		curFace.vecN[0]*curFace.deltaGrad[0]+
		curFace.vecN[1]*curFace.deltaGrad[1]+
		curFace.vecN[2]*curFace.deltaGrad[2];
	      Float GB=
		curFace.vecB[0]*curFace.deltaGrad[0]+
		curFace.vecB[1]*curFace.deltaGrad[1]+
		curFace.vecB[2]*curFace.deltaGrad[2];	      

	      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
		deltaWeight = getPeriodizedDeltaWeight(coord,wf,curFace);
	      
	      if (debug)
		{
		  Float tmp=deltaWeight+(3.0*GB*PB + 2.0*GN*PN + GT*PT)/4.0;
		  std::cout << 
		    "dw =" <<deltaWeight<<"+(3.0*" << GB <<"*"<< PB <<"+"<<
		    "2.0*" << GN << "*" << PN <<"+"<<
		    GT << "*" << PT <<")/4 \n   = "<< tmp<<"\n";
		  
		  Float add=PT*PN*PB*tmp;
		  std::cout << "contrib += " /*<< contrib <<"+" */<<
		    PT<<"*"<<PN<<"*"<<PB<<"*"<<
		    deltaWeight + (3.0*GB*PB + 2.0*GN*PN + GT*PT)/4.0 <<" = "<<
		    contrib <<"+"<< add << " = " << contrib+add <<"\n";
		}
	      
	      deltaWeight += (3.0*GB*PB + 2.0*GN*PN + GT*PT)/4.0;
	    }

	  Float localContrib=PT*PN*PB*deltaWeight;
	  contrib += localContrib;

	  if (fabs(localContrib)>fabs(maxContrib))
	    maxContrib=localContrib;
	  if (fabs(contrib)>fabs(maxContrib))
	    maxContrib=contrib;
	  //contrib +=  PT*PN*PB*deltaWeight;
	  
	  // std::abs not working correctly for Float=float128 ...
	  /*	  
	    {
	      if (localContrib<0)
		{
		  if ((-localContrib)>maxContrib)
		    maxContrib=-localContrib;
		}
	      else if (localContrib>maxContrib)
		maxContrib=localContrib;
	    }
	  */
	}
      //if (debug) std::cout << "Required precision: " << (maxContrib/contrib)  <<".\n";
      return contrib;
    }

    template <class T, class WF, class LFloat=Float>    
    Float computeVoxelFacetSimplexEdgeContrib(const T *coord, 
					      int normalDim, Float normalSign,
					      const WF &wf, const IncidentEdge &edge,
					      Float &maxContrib,bool dbg=false)
    {
      maxContrib=0;
      bool almostDegenerate=(fabs(edge.vecT[normalDim])<1.E-8);
      //if (almostDegenerate) printf("ORIENT almostdegenerate! \n");
      /*
      if (dbg) almostDegenerate=false;
      if (almostDegenerate)
	return hlp::numericStaticCast<Float>
	  (computeVoxelFacetSimplexEdgeContrib<T,WF,HFloat>
	   (coord,normalDim,normalSign,wf,edge,true));
      */
     	
      /*
      if (dbg)
	{
	  printf("Normal[%d] = %Lg -> vecT : %Lg\n",
		 normalDim,nst(normalSign),nst(edge.vecT[normalDim]));
	  printf("Coord = (%20.20Lg %20.20Lg %20.20Lg)\n",
		 nst(coord[0]),nst(coord[1]),nst(coord[2]));
	  if (almostDegenerate) printf("Almost degenerate.\n");
	}
      */

      LFloat contrib=0;      //T coord[NDIM]={coordd[0],coordd[1],coordd[2]};
      // Contribution from the voxel face normal
      LFloat PvB = normalSign*coord[normalDim];
      
      // Part of the contribution from the simplex Edge tangent vector
      LFloat PsT = 
	edge.vecT[0] * coord[0]+
       	edge.vecT[1] * coord[1]+
       	edge.vecT[2] * coord[2]; 
      
      //if (dbg) std::cout << "PST = " <<PsT<<std::endl;
      
      //if ((edge.vecT[normalDim]>0)!=(normalSign>0)) PsT = -PsT;
     
      for (int f=0;f<edge.face.size();++f)
	{ 
	  const IncidentFace &curFace=edge.face[f];
	  //if (curFace.nextWeight == curFace.prevWeight) continue;

	  //if (dbg) std::cout<<"Normal along orth : "<<curFace.vecB[normalDim]<<std::endl;
	  
	  // Part of the contribution from the simplex Edge normal vector
	  LFloat PsN=
	    curFace.vecN[0] * coord[0]+
	    curFace.vecN[1] * coord[1]+
	    curFace.vecN[2] * coord[2];

	  // Part of the contribution from the simplex Edge binormal vector
	  LFloat PsB=
	    curFace.vecB[0] * coord[0]+
	    curFace.vecB[1] * coord[1]+
	    curFace.vecB[2] * coord[2];
	  //if (dbg) std::cout << "PSB/PSN["<<f<<"] = "<<PsN<<" "<<PsB<<std::endl;

	  // Tangent vector along the intersection of the simplex facet and voxel facet	  
	  LFloat vecVT[3];
	  // Normal vector along the voxel facet (projection of simplex's binormal vector)
	  LFloat vecVN[3];
	  // Normal vector along the simplex, when starting along the voxel/simplex facet
	  // intersection
	  //LFloat vecSN[3]; //unused
	  
	  // projection of simplex facet normal onto the voxel face
	  std::copy(curFace.vecB,curFace.vecB+NDIM,vecVN);
	  vecVN[normalDim]=0;
	  meshGeometry->template normalize_noCheck<LFloat,NDIM>(vecVN);
	  
	  LFloat PvT;
	  LFloat PvN;
	  LFloat GvT;
	  LFloat GvN;
	  if (normalDim==0)
	    {	      
	      vecVT[0]=0;
	      vecVT[1]=vecVN[2];
	      vecVT[2]=-vecVN[1];
	      
	      if (almostDegenerate)
		{
		  // This is necessary because vecT may not have enough 
		  // precision to determine the sign !		  
		  LFloat Tsign = meshGeometry->coordDiffSign
		    (edge.vertex->getCoord(normalDim),
		     edge.otherVertex->getCoord(normalDim),
		     normalDim);

		  //(edge.vecT[normalDim]>0)
		  if ((Tsign>0)==(curFace.vecN[normalDim]>0))
		    {
		      if (vecVT[1]*edge.vecT[1]+vecVT[2]*edge.vecT[2]>0)
			{
			  //if (dbg) printf("T was negated1\n");
			  vecVT[1]=-vecVT[1];
			  vecVT[2]=-vecVT[2];
			}
		    }
		  else
		    {
		      if (vecVT[1]*edge.vecT[1]+vecVT[2]*edge.vecT[2]<0)
			{
			  //if (dbg) printf("T was negated2\n");
			  vecVT[1]=-vecVT[1];
			  vecVT[2]=-vecVT[2];
			}
		    }
		}
	      else if (vecVT[1]*curFace.vecN[1]+vecVT[2]*curFace.vecN[2]<0)
		{
		  //if (dbg) printf("T was negated0\n");
		  vecVT[1]=-vecVT[1];
		  vecVT[2]=-vecVT[2];
		}

	      PvT = vecVT[1]*coord[1]+vecVT[2]*coord[2];
	      PvN = vecVN[1]*coord[1]+vecVN[2]*coord[2];	
	      if (WF::ORDER>0)
		{
		  GvT=vecVT[1]*curFace.deltaGrad[1]+vecVT[2]*curFace.deltaGrad[2];
		  GvN=vecVN[1]*curFace.deltaGrad[1]+vecVN[2]*curFace.deltaGrad[2];
		}
	    }
	  else if (normalDim==1)
	    {	      
	      vecVT[0]=-vecVN[2];
	      vecVT[1]=0;
	      vecVT[2]=vecVN[0];
	      
	      if (almostDegenerate)
		{
		  if ((edge.vecT[normalDim]>0)==(curFace.vecN[normalDim]>0))
		    {
		      if (vecVT[0]*edge.vecT[0]+vecVT[2]*edge.vecT[2]>0)
			{
			  vecVT[0]=-vecVT[0];
			  vecVT[2]=-vecVT[2];
			}
		    }
		  else
		    {
		      if (vecVT[0]*edge.vecT[0]+vecVT[2]*edge.vecT[2]<0)
			{
			  vecVT[0]=-vecVT[0];
			  vecVT[2]=-vecVT[2];
			}
		    }
		}
	      else if (vecVT[0]*curFace.vecN[0]+vecVT[2]*curFace.vecN[2]<0)
		{
		  vecVT[0]=-vecVT[0];
		  vecVT[2]=-vecVT[2];
		}

	      PvT = vecVT[0]*coord[0]+vecVT[2]*coord[2];
	      PvN = vecVN[0]*coord[0]+vecVN[2]*coord[2];	
	      if (WF::ORDER>0)
		{
		  GvT=vecVT[0]*curFace.deltaGrad[0]+vecVT[2]*curFace.deltaGrad[2];
		  GvN=vecVN[0]*curFace.deltaGrad[0]+vecVN[2]*curFace.deltaGrad[2];
		}
	    }
	  else
	    {
	      vecVT[0]=vecVN[1];
	      vecVT[1]=-vecVN[0];
	      vecVT[2]=0;

	      if (almostDegenerate)
		{
		  if ((edge.vecT[normalDim]>0)==(curFace.vecN[normalDim]>0))
		    {
		      if (vecVT[0]*edge.vecT[0]+vecVT[1]*edge.vecT[1]>0)
			{
			  vecVT[0]=-vecVT[0];
			  vecVT[1]=-vecVT[1];
			}
		    }
		  else
		    {
		      if (vecVT[0]*edge.vecT[0]+vecVT[1]*edge.vecT[1]<0)
			{
			  vecVT[0]=-vecVT[0];
			  vecVT[1]=-vecVT[1];
			}
		    }
		}
	      else if (vecVT[0]*curFace.vecN[0]+vecVT[1]*curFace.vecN[1]<0)
		{
		  vecVT[0]=-vecVT[0];
		  vecVT[1]=-vecVT[1];
		}

	      PvT = vecVT[0]*coord[0]+vecVT[1]*coord[1];
	      PvN = vecVN[0]*coord[0]+vecVN[1]*coord[1];
	      if (WF::ORDER>0)
		{
		  GvT=vecVT[0]*curFace.deltaGrad[0]+vecVT[1]*curFace.deltaGrad[1];
		  GvN=vecVN[0]*curFace.deltaGrad[0]+vecVN[1]*curFace.deltaGrad[1];
		}
	    }

	  // vecSVN is orthogonal to the tangent along the simplex/voxel facets and to 
	  // the simplex facet normal (i.e. binormal vector vecB)
	  LFloat vecSVN[3];
	  vecSVN[0]=vecVT[1]*curFace.vecB[2]-vecVT[2]*curFace.vecB[1];
	  vecSVN[1]=vecVT[2]*curFace.vecB[0]-vecVT[0]*curFace.vecB[2];
	  vecSVN[2]=vecVT[0]*curFace.vecB[1]-vecVT[1]*curFace.vecB[0];
	  
	  LFloat PsvN =
	    vecSVN[0]*coord[0]+
	    vecSVN[1]*coord[1]+
	    vecSVN[2]*coord[2];

	  //if (dbg) std::cout << "PvT/PvN  :: GvT/GvN = " << PvT<<" "<<PvN<<" :: "<<GvT<<" "<<GvN<<std::endl;
	  //if (dbg) std::cout << "PsvN = " << PsvN<<" ("<<vecSVN[0]<<" "<<vecSVN[1]<<" "<<vecSVN[2]<<")"<<std::endl;
	  //if ((vecSVN[normalDim]>0)!=(normalSign>0)) PsvN = -PsvN;

	  // Contrib from tangent along voxel / normal along simplex / binormal of simplex
	  //LFloat localContrib = vTContrib*svNContrib*sBContrib;	  

	  // Contrib from tangent along voxel / normal along voxel / 
	  // binormal orthogonal to voxel face
	  //localContrib += vTContrib*vNContrib*vBContrib;	  

	  // and contrib along simplex's tangent/Normal/binormal
	  //localContrib += sTContrib*sNContrib*sBContrib;
	  Float thisContrib;
	  if (WF::ORDER>0)
	    {
	      LFloat GsvN=
		vecSVN[0]*curFace.deltaGrad[0]+
		vecSVN[1]*curFace.deltaGrad[1]+
		vecSVN[2]*curFace.deltaGrad[2];	 
	      LFloat GvB=normalSign*curFace.deltaGrad[normalDim];
	      LFloat GsB=
		curFace.vecB[0] * curFace.deltaGrad[0]+
		curFace.vecB[1] * curFace.deltaGrad[1]+
		curFace.vecB[2] * curFace.deltaGrad[2];
	      LFloat GsN=
		curFace.vecN[0] * curFace.deltaGrad[0]+
		curFace.vecN[1] * curFace.deltaGrad[1]+
		curFace.vecN[2] * curFace.deltaGrad[2];
	      LFloat GsT=
		edge.vecT[0] * curFace.deltaGrad[0]+
		edge.vecT[1] * curFace.deltaGrad[1]+
		edge.vecT[2] * curFace.deltaGrad[2];

	      //if (dbg) std::cout<<"GsvN GvB GsB GsN GsT = "<< GsvN<<" "<<GvB<<" "<<GsB<<" "<<GsN<<" "<<GsT<<std::endl;
		

	      if ((edge.vecT[normalDim]>0)!=(normalSign>0))
		{
		  GsT = -GsT;
		  PsT = -PsT;
		  //if (dbg) printf("PsT  was negated\n");
		}
	      if ((vecSVN[normalDim]>0)!=(normalSign>0)) 
		{
		  PsvN = -PsvN;
		  GsvN = -GsvN;
		  //if (dbg) printf("PsvN  was negated\n");
		}
	      
	      LFloat deltaWeight = curFace.deltaWeight;
	      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
		deltaWeight = getPeriodizedDeltaWeight(coord,wf,curFace);

	      if (fabs(vecSVN[normalDim])<1.E-15)
		{
		  //printf("orient: vecSVN degene !\n");
		  //if (dbg) std::cout <<"Degenerate CASE <1.E-15"<< std::endl;

		  // This may happen when edge.vecT is almost orthogonal to the voxel 
		  // facet normal. Then in some casse, simplex facets may become tangent to
		  // voxel facets, and their intersection cannot be computed safely ... 
		  // Fortunately, in such case, only the simplex edge T/N/B contributes.
		  //if (dbg)
		  //LFloat before = contrib;
		  thisContrib=PsT*PsN*PsB*
		    (deltaWeight + (3.0*GsB*PsB + 2.0*GsN*PsN + GsT*PsT)/4.0);
		  contrib += thisContrib;
		  if (fabs(thisContrib)>fabs(maxContrib))
		    maxContrib=thisContrib;
		  if (fabs(contrib)>fabs(maxContrib))
		    maxContrib=contrib;
		  /*
		  if (dbg)
		    {
		      std::cout << "contrib+=" << PsT*PsN*PsB*
		    (deltaWeight + (3.0*GsB*PsB + 2.0*GsN*PsN + GsT*PsT)/4.0)<< " -> c=" 
				<< contrib <<std::endl;
		    }
		  */
		  
		}
	      else
		{
		  /*
		  if (dbg)
		    {
		      std::cout <<"Regular contrib: " 
				<< PvT*PsvN*PsB*(deltaWeight + (3.0*GsB*PsB + 2.0*GsvN*PsvN + GvT*PvT)/4.0) 
				<< std::endl;
		      std::cout <<"Additional contrib: "
				<< PvT*PvN*PvB*(deltaWeight + (3.0*GvB*PvB + 2.0*GvN*PvN + GvT*PvT)/4.0)+PsT*PsN*PsB*(deltaWeight + (3.0*GsB*PsB + 2.0*GsN*PsN + GsT*PsT)/4.0) 
				<< std::endl;
		    }
		  */
		  
		  //LFloat before = contrib;
		  thisContrib=PvT*PsvN*PsB*
		    (deltaWeight + (3.0*GsB*PsB + 2.0*GsvN*PsvN + GvT*PvT)/4.0);
		  contrib += thisContrib;
		  if (fabs(thisContrib)>fabs(maxContrib))
		    maxContrib=thisContrib;
		  if (fabs(contrib)>fabs(maxContrib))
		    maxContrib=contrib;

		  thisContrib=PvT*PvN*PvB*
		    (deltaWeight + (3.0*GvB*PvB + 2.0*GvN*PvN + GvT*PvT)/4.0);
		  contrib += thisContrib;
		  if (fabs(thisContrib)>fabs(maxContrib))
		    maxContrib=thisContrib;
		  if (fabs(contrib)>fabs(maxContrib))
		    maxContrib=contrib;

		  thisContrib = PsT*PsN*PsB*
		    (deltaWeight + (3.0*GsB*PsB + 2.0*GsN*PsN + GsT*PsT)/4.0);
		  contrib+=thisContrib;
		  if (fabs(thisContrib)>fabs(maxContrib))
		    maxContrib=thisContrib;
		  if (fabs(contrib)>fabs(maxContrib))
		    maxContrib=contrib;
		 
		  //if (dbg) std::cout << " -> c=" << contrib << std::endl;
		  
		}
	    }
	  else
	    {
	      if ((edge.vecT[normalDim]>0)!=(normalSign>0))
		PsT = -PsT;
	      if ((vecSVN[normalDim]>0)!=(normalSign>0)) 	      
		PsvN = -PsvN;
	      
	      if (fabs(vecSVN[normalDim])<1.E-15)/*1E-15)*/
		{	
		  // This may happen when edge.vecT is almost orthogonal to the voxel 
		  // facet normal. Then in some cases, simplex facets may become tangent to
		  // voxel facets, and their intersection cannot be computed safely ... 
		  // Fortunately, in such case, only the simplex edge T/N/B contributes.
		  contrib += curFace.deltaWeight * (PsT*PsN*PsB);		 
		}
	      else
		{
		  thisContrib= curFace.deltaWeight * 
		    (PvT*PsvN*PsB + PvT*PvN*PvB + PsT*PsN*PsB);
		  contrib += thisContrib;
		  if (fabs(thisContrib)>fabs(maxContrib))
		    maxContrib=thisContrib;
		  if (fabs(contrib)>fabs(maxContrib))
		    maxContrib=contrib;
		}
	      /*
	       if (dbg)
		    {
		   printf("(edge.vecT[%d]=%Lg>0)= (%d == %d)\n",normalDim,
			  nst(edge.vecT[normalDim]),(edge.vecT[normalDim]>0),
			  normalSign>0);
		   printf("(vecSVN[%d]=%Lg>0)= (%d == %d)\n",normalDim,
			  nst(vecSVN[normalDim]),(vecSVN[normalDim]>0),
			  normalSign>0);
		   printf("vecSVN=(%Lg,%Lg,%Lg)\n",
			  nst(vecSVN[0]),nst(vecSVN[1]),nst(vecSVN[2]));
		  
		   printf("(%Lg*%Lg*%Lg=%Lg) + (%Lg*%Lg*%Lg=%Lg) + (%Lg*%Lg*%Lg=%Lg) = %Lg\n",
			  nst(PvT),nst(PsvN),nst(PsB),nst(PvT*PsvN*PsB), 
			  nst(PvT),nst(PvN),nst(PvB),nst(PvT*PvN*PvB), 
			  nst(PsT),nst(PsN),nst(PsB),nst(PsT*PsN*PsB),
			  nst(PvT*PsvN*PsB + PvT*PvN*PvB + PsT*PsN*PsB));
		    }
	      */
	    }
	}

      return hlp::numericStaticCast<Float>(contrib);
    }

    
    // projectedFacetNormal is the normalized normal to the facet. Its direction
    // is not relevant (it can be pointing inward or outward simplex)

    // NB: the voxel edge orientation (in vEdgeContribSign[dim]) is modified
    // by this function so that it becomes aligned with the facetNormal
    // NB2: the order of segNeighbors must be such that segNeighbor[0] and
    // segNeighbor[3] do not share a face
    template <class WF>
    Float computeVoxelEdgeSimplexFacetContrib(Simplex *simplex, 
					      int simplexNeighborIndex, 
					      Float projectedFacetNormal[NDIM], 
					      int dim,
					      Float vEdgeContribSign[NDIM], 
					      Float intersection[NDIM],
					      const WF &wf, 
					      Float &maxContrib,
					      bool coordsAreConsistent=false,
					      bool check=false) const
    { 
      maxContrib=0;
      if (check) std::cout<<"Facet normal : "<<projectedFacetNormal[dim]<<std::endl;

      //if (fabs(projectedFacetNormal[dim])<1.E-4) return 0;

      //int nProcessed=0;
      Simplex *nei = simplex->getNeighbor(simplexNeighborIndex);
      
      //Float wOut=validSimplex(nei)?wf.get(nei):0;      
      //Float wIn=validSimplex(simplex)?wf.get(simplex):0;      

      Float wIn;
      Float wOut;
      Float gradIn[NDIM];
      Float gradOut[NDIM];
      
      // Compute the weights
      getPeriodizedWeightAndGrad(intersection,wOut,gradOut,wf,nei);
      getPeriodizedWeightAndGrad(intersection,wIn,gradIn,wf,simplex);
      if (check) 
	{
	  std::cout << "wei in/out:"<<wIn<<" "<<wOut<<std::endl;
	  std::cout << "grad in : "<<gradIn[0]<<" "<<gradIn[1]<<" "<<gradIn[2]<<std::endl;
	  std::cout << "grad out: "<<gradOut[0]<<" "<<gradOut[1]<<" "<<gradOut[2]<<std::endl;
	}
      // First check the orientation of the normal
      // => make it point OUTWARD
      // FIXME: What about when the simplex is degenerate ? then we are screwed ...
      Float tmpVec[NDIM];
      Float vc[NDIM];
      for (int i=0;i<NDIM;++i) vc[i]=simplex->getVertex(simplexNeighborIndex)->getCoord(i);
      meshGeometry->template getVector<Float,Float,NDIM>(intersection,vc,tmpVec);
		
      Float dot = meshGeometry->template 
	dot_noCheck<Float,Float,NDIM>(tmpVec,projectedFacetNormal);               

      Float facetNormal[NDIM];
      if (check) std::cout << "dot1="<<dot<<std::endl;
      if (dot>0) 
	{
	  // Normal has the wrong orientation
	  facetNormal[0]=-projectedFacetNormal[0];
	  facetNormal[1]=-projectedFacetNormal[1];
	  facetNormal[2]=-projectedFacetNormal[2];
	}
      else std::copy(projectedFacetNormal,projectedFacetNormal+NDIM,facetNormal);

      // and check whether we are folding along this facet     
      if (validSimplex(nei))
	{
	  Vertex *oppVertex=nei->getVertex(0);
	  for (int i=1;i<Simplex::NNEI;++i)
	    {
	      if (nei->getNeighbor(i)==simplex)
		oppVertex = nei->getVertex(i);
	    }
	  for (int i=0;i<NDIM;++i) vc[i]=oppVertex->getCoord(i);
	  meshGeometry->template 
	    getVector<Float,Float,NDIM>(intersection,vc,tmpVec);
	  
	  Float dot = meshGeometry->template 
	    dot_noCheck<Float,Float,NDIM>(tmpVec,facetNormal);

	  // We are folding => neighbor contributes with the same sign as simplex	  
	  if (dot<0) 
	    {
	      wOut=-wOut;
	      if (WF::ORDER>0)
		{
		  gradOut[0]=-gradOut[0];
		  gradOut[1]=-gradOut[1];
		  gradOut[2]=-gradOut[2];
		}
	    }
	}

      Float deltaGrad[NDIM];
      if (WF::ORDER>0)
	{
	  for (int i=0;i<NDIM;++i) 
	    deltaGrad[i]=gradOut[i]-gradIn[i];
	}

      //if (wOut == wIn) return 0;
      
      // Contribution from the normal of the simplex facet
      Float PsB=meshGeometry->template dot_noCheck<Float,Float,NDIM>
	(facetNormal,intersection);
      Float GsB;
      if (WF::ORDER>0)
	{
	  GsB = 
	    facetNormal[0]*deltaGrad[0]+
	    facetNormal[1]*deltaGrad[1]+
	    facetNormal[2]*deltaGrad[2];
	}      
         	  
      // The 2 tangent vectors at the simplex facet / voxel facet intersection
      Float vecVT[2][3];
      // The 2 normal vectors of vecVT along the voxel facets
      Float vecVN[2][3];
      // The 2 normal vectors of vecVT along the simplex facet
      Float vecSN[2][3];

      // Contribution from voxel facets normals
      Float PvB[2];
      // Contribution from the vecVT
      Float PvT[2];	  
      // Contribution from vecVN
      Float PvN[2];
      // Contribution from vecSN
      Float PsN[2];

      Float GvB[2];
      Float GvT[2];
      Float GvN[2];
    
      // We need to have the contribution along dim in the same direction
      // for the facet normal and the voxel edge      
      if (facetNormal[dim]>0) 
	vEdgeContribSign[dim]=1; 
      else 
	vEdgeContribSign[dim]=-1;

      // each vecVN is computed from the facet normal projected onto planes orthogonal
      // to 'dim'.
      std::copy(facetNormal,facetNormal+NDIM,vecVN[0]);
      std::copy(facetNormal,facetNormal+NDIM,vecVN[1]);

      if (dim==0)
	{	   
	  vecVN[0][1]=0; // orthogonal to Y
	  vecVN[1][2]=0; // orthogonal to Z
	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[0]);	    
	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[1]);

	  // vector orthogonal to simplex and voxel facets normals
	  vecVT[0][0]=-vecVN[0][2];
	  vecVT[0][1]=0;
	  vecVT[0][2]=vecVN[0][0];
	  vecVT[1][0]=vecVN[1][1];
	  vecVT[1][1]=-vecVN[1][0];
	  vecVT[1][2]=0;
	  
	  if (check) std::cout<<"vecVT[0] = "<<vecVT[0][2]<<" "<<vecVT[1][1]<<std::endl;
	  //if (check) printf("vecVT[0] = %e %e\n",vecVT[0][2],vecVT[1][1]);
	  // Correct orientation if needed
	  bool swap=false;
	  if (fabs(vecVT[0][2])<1.E-10)
	    {
	      //printf("ORIENT02! \n");
	      // In that case we lack precision to orient VT reliably using the 
	      // facet normal component along 2, but we can use its component along dim ...
	      if ((facetNormal[2]>0)==(vEdgeContribSign[2]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[0][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[0][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[2]>0)!=(vecVT[0][2]>0)) swap=true;

	  if (swap)
	    {
	      vecVT[0][0]=-vecVT[0][0];
	      vecVT[0][2]=-vecVT[0][2];
	    }
	  
	  swap=false;
	  if (fabs(vecVT[1][1])<1.E-10)
	    {	
	      //printf("ORIENT11! \n");
	      // In that case we lack precision to orient VT reliably using the 
	      // facet normal component along 1, but we can use its component along dim ...
	      if ((facetNormal[1]>0)==(vEdgeContribSign[1]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[1][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[1][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[1]>0)!=(vecVT[1][1]>0)) swap=true;

	  if (swap)
	    {
	      vecVT[1][0]=-vecVT[1][0];
	      vecVT[1][1]=-vecVT[1][1];
	    }
					    
	  // and compute respective contribution
	  PvB[0]=vEdgeContribSign[1]*intersection[1];
	  PvB[1]=vEdgeContribSign[2]*intersection[2];

	  PvT[0]=
	    intersection[0]*vecVT[0][0]+
	    intersection[2]*vecVT[0][2];
	  PvT[1]=
	    intersection[0]*vecVT[1][0]+
	    intersection[1]*vecVT[1][1];

	  PvN[0]=
	    intersection[0]*vecVN[0][0]+
	    intersection[2]*vecVN[0][2];
	  PvN[1]=
	    intersection[0]*vecVN[1][0]+
	    intersection[1]*vecVN[1][1];

	  if (WF::ORDER>0)
	    {
	      GvB[0]=vEdgeContribSign[1]*deltaGrad[1];
	      GvB[1]=vEdgeContribSign[2]*deltaGrad[2];
	      GvT[0]=
		deltaGrad[0]*vecVT[0][0]+
		deltaGrad[2]*vecVT[0][2];
	      GvT[1]=
		deltaGrad[0]*vecVT[1][0]+
		deltaGrad[1]*vecVT[1][1];
	      GvN[0]=
		deltaGrad[0]*vecVN[0][0]+
		deltaGrad[2]*vecVN[0][2];
	      GvN[1]=
		deltaGrad[0]*vecVN[1][0]+
		deltaGrad[1]*vecVN[1][1];
	    }
					    
	  vecSN[0][0]=vecVT[0][1]*facetNormal[2]-vecVT[0][2]*facetNormal[1];
	  vecSN[0][1]=vecVT[0][2]*facetNormal[0]-vecVT[0][0]*facetNormal[2];
	  vecSN[0][2]=vecVT[0][0]*facetNormal[1]-vecVT[0][1]*facetNormal[0];
	  PsN[0] = 
	    vecSN[0][0]*intersection[0]+
	    vecSN[0][1]*intersection[1]+
	    vecSN[0][2]*intersection[2];
	  if ((vEdgeContribSign[1]>0)!=(vecSN[0][1]>0))
	    {
	      PsN[0]=-PsN[0];
	      for (int k=0;k<NDIM;++k) vecSN[0][k]=-vecSN[0][k];
	    }
					    
	  vecSN[1][0]=vecVT[1][1]*facetNormal[2]-vecVT[1][2]*facetNormal[1];
	  vecSN[1][1]=vecVT[1][2]*facetNormal[0]-vecVT[1][0]*facetNormal[2];
	  vecSN[1][2]=vecVT[1][0]*facetNormal[1]-vecVT[1][1]*facetNormal[0];
	  PsN[1] = 
	    vecSN[1][0]*intersection[0]+
	    vecSN[1][1]*intersection[1]+
	    vecSN[1][2]*intersection[2];
	  if ((vEdgeContribSign[2]>0)!=(vecSN[1][2]>0))
	    {
	      PsN[1]=-PsN[1];
	      for (int k=0;k<NDIM;++k) vecSN[1][k]=-vecSN[1][k];
	    }
	  if (check) std::cout<<"vecSN[0] = "<<vecSN[0][1]<<" "<<vecSN[1][2]<<std::endl;
	  //if (check) printf("vecSN[0] = %e %e\n",vecSN[0][1],vecSN[1][2]);
	}
      else if (dim==1)
	{
	  vecVN[0][0]=0;
	  vecVN[1][2]=0;

	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[0]);	    
	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[1]);

	  vecVT[0][0]=0;
	  vecVT[0][1]=vecVN[0][2];
	  vecVT[0][2]=-vecVN[0][1];
	  vecVT[1][0]=vecVN[1][1];
	  vecVT[1][1]=-vecVN[1][0];
	  vecVT[1][2]=0;

	  if (check) std::cout<<"vecVT[1] = "<<vecVT[0][2]<<" "<<vecVT[1][0]<<std::endl;
	  //if (check) printf("vecVT[1] = %e %e\n",vecVT[0][2],vecVT[1][0]);
	  
	  bool swap=false;
	  if (fabs(vecVT[0][2])<1.E-10)
	    {
	      //printf("1ORIENT02! \n");
	      if ((facetNormal[2]>0)==(vEdgeContribSign[2]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[0][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[0][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[2]>0)!=(vecVT[0][2]>0)) swap=true;
	  
	  if (swap)
	    {
	      vecVT[0][1]=-vecVT[0][1];
	      vecVT[0][2]=-vecVT[0][2];
	    }

	  
	  swap=false;
	  if (fabs(vecVT[1][0])<1.E-10)
	    {
	      //printf("1ORIENT10! \n");
	      if ((facetNormal[0]>0)==(vEdgeContribSign[0]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[1][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[1][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[0]>0)!=(vecVT[1][0]>0)) swap=true;

	  if (swap)
	    {
	      vecVT[1][0]=-vecVT[1][0];
	      vecVT[1][1]=-vecVT[1][1];
	    }

	  PvB[0]=vEdgeContribSign[0]*intersection[0];
	  PvB[1]=vEdgeContribSign[2]*intersection[2];
					    
	  PvT[0]=
	    intersection[1]*vecVT[0][1]+
	    intersection[2]*vecVT[0][2];
	  PvT[1]=
	    intersection[0]*vecVT[1][0]+
	    intersection[1]*vecVT[1][1];

	  PvN[0]=
	    intersection[1]*vecVN[0][1]+
	    intersection[2]*vecVN[0][2];
	  PvN[1]=
	    intersection[0]*vecVN[1][0]+
	    intersection[1]*vecVN[1][1];

	  if (WF::ORDER>0)
	    {
	      GvB[0]=vEdgeContribSign[0]*deltaGrad[0];
	      GvB[1]=vEdgeContribSign[2]*deltaGrad[2];
	      GvT[0]=
		deltaGrad[1]*vecVT[0][1]+
		deltaGrad[2]*vecVT[0][2];
	      GvT[1]=
		deltaGrad[0]*vecVT[1][0]+
		deltaGrad[1]*vecVT[1][1];
	      GvN[0]=
		deltaGrad[1]*vecVN[0][1]+
		deltaGrad[2]*vecVN[0][2];
	      GvN[1]=
		deltaGrad[0]*vecVN[1][0]+
		deltaGrad[1]*vecVN[1][1];
	    }

	  vecSN[0][0]=vecVT[0][1]*facetNormal[2]-vecVT[0][2]*facetNormal[1];
	  vecSN[0][1]=vecVT[0][2]*facetNormal[0]-vecVT[0][0]*facetNormal[2];
	  vecSN[0][2]=vecVT[0][0]*facetNormal[1]-vecVT[0][1]*facetNormal[0];
	  PsN[0] = 
	    vecSN[0][0]*intersection[0]+
	    vecSN[0][1]*intersection[1]+
	    vecSN[0][2]*intersection[2];
	  if ((vEdgeContribSign[0]>0)!=(vecSN[0][0]>0))
	    {
	      PsN[0]=-PsN[0];
	      for (int k=0;k<NDIM;++k) vecSN[0][k]=-vecSN[0][k];
	    }

	  vecSN[1][0]=vecVT[1][1]*facetNormal[2]-vecVT[1][2]*facetNormal[1];
	  vecSN[1][1]=vecVT[1][2]*facetNormal[0]-vecVT[1][0]*facetNormal[2];
	  vecSN[1][2]=vecVT[1][0]*facetNormal[1]-vecVT[1][1]*facetNormal[0];
	  PsN[1] = 
	    vecSN[1][0]*intersection[0]+
	    vecSN[1][1]*intersection[1]+
	    vecSN[1][2]*intersection[2];
	  if ((vEdgeContribSign[2]>0)!=(vecSN[1][2]>0))
	    {
	      PsN[1]=-PsN[1];
	      for (int k=0;k<NDIM;++k) vecSN[1][k]=-vecSN[1][k];
	    }	 

	  if (check) std::cout<<"vecSN[1]"<<vecSN[0][0]<<" "<<vecSN[1][2]<<std::endl;
	  //if (check) printf("vecSN[1] = %e %e\n",vecSN[0][0],vecSN[1][2]);
	}
      else
	{
	  vecVN[0][0]=0; // orthogonal to X
	  vecVN[1][1]=0; // orthogonal to Y
	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[0]);
	  meshGeometry->template normalize_noCheck<Float,NDIM>(vecVN[1]);

	  vecVT[0][0]=0;
	  vecVT[0][1]=vecVN[0][2];
	  vecVT[0][2]=-vecVN[0][1];
	  vecVT[1][0]=-vecVN[1][2];
	  vecVT[1][1]=0;
	  vecVT[1][2]=vecVN[1][0];

	  if (check) std::cout<<"vecVT[2]"<<vecVT[0][1]<<" "<<vecVT[1][0]<<std::endl;
	  //if (check) printf("vecVT[2] = %e %e\n",vecVT[0][1],vecVT[1][0]);
	  
	  bool swap=false;
	  if (fabs(vecVT[0][1])<1.E-10)
	    {	      
	      //printf("2ORIENT01! \n");
	      if ((facetNormal[1]>0)==(vEdgeContribSign[1]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[0][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[0][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[1]>0)!=(vecVT[0][1]>0)) swap=true;

	  if (swap)
	    {
	      vecVT[0][1]=-vecVT[0][1];
	      vecVT[0][2]=-vecVT[0][2];
	    }
	  
	  swap=false;
	  if (fabs(vecVT[1][0])<1.E-10)
	    {
	      //printf("2ORIENT10! \n");
	      //if (vEdgeContribSign[dim]>0)
	      if ((facetNormal[0]>0)==(vEdgeContribSign[0]>0))
		{
		  if ((vEdgeContribSign[dim]>0)==(vecVT[1][dim]>0))
		    swap=true;
		}
	      else
		{
		  if ((vEdgeContribSign[dim]>0)!=(vecVT[1][dim]>0))
		    swap=true;
		}
	    }
	  else if ((vEdgeContribSign[0]>0)!=(vecVT[1][0]>0)) swap=true;

	  if (swap)
	    {
	      vecVT[1][0]=-vecVT[1][0];
	      vecVT[1][2]=-vecVT[1][2];
	    }

	  PvB[0]=vEdgeContribSign[0]*intersection[0];
	  PvB[1]=vEdgeContribSign[1]*intersection[1];

	  PvT[0]=
	    intersection[1]*vecVT[0][1]+
	    intersection[2]*vecVT[0][2];
	  PvT[1]=
	    intersection[0]*vecVT[1][0]+
	    intersection[2]*vecVT[1][2];

	  PvN[0]=
	    intersection[1]*vecVN[0][1]+
	    intersection[2]*vecVN[0][2];
	  PvN[1]=
	    intersection[0]*vecVN[1][0]+
	    intersection[2]*vecVN[1][2];

	  if (WF::ORDER>0)
	    {
	      GvB[0]=vEdgeContribSign[0]*deltaGrad[0];
	      GvB[1]=vEdgeContribSign[1]*deltaGrad[1];
	      GvT[0]=
		deltaGrad[1]*vecVT[0][1]+
		deltaGrad[2]*vecVT[0][2];
	      GvT[1]=
		deltaGrad[0]*vecVT[1][0]+
		deltaGrad[2]*vecVT[1][2];
	      GvN[0]=
		deltaGrad[1]*vecVN[0][1]+
		deltaGrad[2]*vecVN[0][2];
	      GvN[1]=
		deltaGrad[0]*vecVN[1][0]+
		deltaGrad[2]*vecVN[1][2];
	    }

	  vecSN[0][0]=vecVT[0][1]*facetNormal[2]-vecVT[0][2]*facetNormal[1];
	  vecSN[0][1]=vecVT[0][2]*facetNormal[0]-vecVT[0][0]*facetNormal[2];
	  vecSN[0][2]=vecVT[0][0]*facetNormal[1]-vecVT[0][1]*facetNormal[0];
	  PsN[0] = 
	    vecSN[0][0]*intersection[0]+
	    vecSN[0][1]*intersection[1]+
	    vecSN[0][2]*intersection[2];
	  if ((vEdgeContribSign[0]>0)!=(vecSN[0][0]>0))
	    {
	      PsN[0]=-PsN[0];
	      for (int k=0;k<NDIM;++k) vecSN[0][k]=-vecSN[0][k];
	    }

	  vecSN[1][0]=vecVT[1][1]*facetNormal[2]-vecVT[1][2]*facetNormal[1];
	  vecSN[1][1]=vecVT[1][2]*facetNormal[0]-vecVT[1][0]*facetNormal[2];
	  vecSN[1][2]=vecVT[1][0]*facetNormal[1]-vecVT[1][1]*facetNormal[0];
	  PsN[1] = 
	    vecSN[1][0]*intersection[0]+
	    vecSN[1][1]*intersection[1]+
	    vecSN[1][2]*intersection[2];
	  if ((vEdgeContribSign[1]>0)!=(vecSN[1][1]>0))
	    {
	      PsN[1]=-PsN[1];
	      for (int k=0;k<NDIM;++k) vecSN[1][k]=-vecSN[1][k];
	    }
	  if (check) std::cout<<"vecSN[2]"<<vecSN[0][0]<<" "<<vecSN[1][1]<<std::endl;
	  //if (check) printf("vecSN[2] = %e %e\n",vecSN[0][0],vecSN[1][1]);
	}      

      if (check) 
	{
	  std::cout << "Contrib = (" 
		    << PvT[0]<<" "<<PsN[0]<<" "
		    << PvT[1]<<" "<<PsN[1]<<" "<<PsB<<")"
		    <<std::endl;
	  // printf("contrib = (%Lg %Lg %Lg %Lg +%Lg)\n",
	  // 		PvT[0],PsN[0],PvT[1],PsN[1],PsB);
	}
      
      Float contrib;
      if (WF::ORDER>0)
	{
	  Float deltaWeight = wOut-wIn;
	  Float GveT=vEdgeContribSign[dim]*deltaGrad[dim];
	  Float PveT=vEdgeContribSign[dim]*intersection[dim];
	  Float edgeProduct=vEdgeContribSign[0]*vEdgeContribSign[1]*vEdgeContribSign[2];
	  contrib = edgeProduct*intersection[0]*intersection[1]*intersection[2];	  

	  contrib*=2.0*deltaWeight + 
	    (5.0*(vEdgeContribSign[0]*intersection[0]*
		   vEdgeContribSign[0]*deltaGrad[0] +
		   vEdgeContribSign[1]*intersection[1]*
		   vEdgeContribSign[1]*deltaGrad[1] +
		   vEdgeContribSign[2]*intersection[2]*
		   vEdgeContribSign[2]*deltaGrad[2]) -
	     3.0*GveT*PveT)/4.0;
	  if (fabs(contrib)>fabs(maxContrib)) maxContrib=contrib;

	  Float GsN[2];
	  GsN[0]=
	    vecSN[0][0]*deltaGrad[0]+
	    vecSN[0][1]*deltaGrad[1]+
	    vecSN[0][2]*deltaGrad[2];
	  GsN[1]=
	    vecSN[1][0]*deltaGrad[0]+
	    vecSN[1][1]*deltaGrad[1]+
	    vecSN[1][2]*deltaGrad[2];

	  Float thisContrib=PvT[0]*PsN[0]*PsB* 
	    (deltaWeight+(3.0*GsB*PsB + 2.0*GsN[0]*PsN[0] + GvT[0]*PvT[0])/4.0);
	  contrib += thisContrib;
	  if (fabs(thisContrib)>fabs(maxContrib)) maxContrib=contrib;
	  if (fabs(contrib)>fabs(maxContrib)) maxContrib=contrib;
	    
	  thisContrib=PvT[1]*PsN[1]*PsB* 
	    (deltaWeight+(3.0*GsB*PsB + 2.0*GsN[1]*PsN[1] + GvT[1]*PvT[1])/4.0);
	  contrib += thisContrib;
	  if (fabs(thisContrib)>fabs(maxContrib)) maxContrib=contrib;
	  if (fabs(contrib)>fabs(maxContrib)) maxContrib=contrib;

	  thisContrib=PvT[0]*PvN[0]*PvB[0]*
	    (deltaWeight+(3.0*GvB[0]*PvB[0] + 2.0*GvN[0]*PvN[0] + GvT[0]*PvT[0])/4.0);
	  contrib += thisContrib;
	  if (fabs(thisContrib)>fabs(maxContrib)) maxContrib=contrib;
	  if (fabs(contrib)>fabs(maxContrib)) maxContrib=contrib;

	  thisContrib=PvT[1]*PvN[1]*PvB[1]*
	    (deltaWeight+(3.0*GvB[1]*PvB[1] + 2.0*GvN[1]*PvN[1] + GvT[1]*PvT[1])/4.0);
	  contrib += thisContrib;
	  if (fabs(thisContrib)>fabs(maxContrib)) maxContrib=contrib;
	  if (fabs(contrib)>fabs(maxContrib)) maxContrib=contrib;
	}
      else
	{	  
	  // First add the contrib that starts along the voxel edge,
	  // there are 2 of them !					
	  contrib= 2.0*
	    vEdgeContribSign[0]*vEdgeContribSign[1]*vEdgeContribSign[2]*
	    intersection[0]*intersection[1]*intersection[2];      
	  
	  // Add the contribution from the simplex/voxel face intersection
	  // along the simplex face
	  contrib += 
	    PvT[0]*PsN[0]*PsB+
	    PvT[1]*PsN[1]*PsB;

	  // and compute the contrib from the simplex/voxel face intersection
	  // along the voxel face					
	  contrib += 
	    PvT[0]*PvB[0]*PvN[0]+
	    PvT[1]*PvB[1]*PvN[1];
	  
	  // Finaly get the total contribution for the two sides of the simplex
	  contrib *= wOut-wIn;
	}

      //if (contrib!=contrib) throw 0;

      /*
      nProcessed+=4;
      (*out)=std::make_pair(segNeighbor[0],contrib);++out;
      (*out)=std::make_pair(segNeighbor[1],-contrib);++out;
      (*out)=std::make_pair(segNeighbor[2],-contrib);++out;	  
      (*out)=std::make_pair(segNeighbor[3],contrib);++out;
      */
#ifdef DEBUGDUMP
      if (findVertex(0.248461,0.0859375,0.171875,intersection))
	{
	  printf("Intersection @(%g %g %g) -> C= %g (w=%e)\n",
		 intersection[0],intersection[1],intersection[2],contrib,
		 (wOut-wIn));	    

	  double fac = 1.0/128.0;
	  Coord c[Simplex::NDIM_W];
	  fprintf(fl[0],"ANDNET\n3\n%d\n",1*4+2*7);
	  //for (int k=0;k<Simplex::Facet::NVERT;++k)
	  //PRINTVERTEX(simplex->getFacetHandle(f)->getVertex(k));
	  for (int k=0;k<Simplex::NVERT;++k)
	    PRINTVERTEX(simplex->getVertex(k));
					
					
	  double cs[2][3]={{0,0,0},{0,0,0}};
	  for (int j=0;j<2;++j)
	    {
	      PRINTCOORD(intersection);
	      for (int k=0;k<NDIM;++k) c[k]=intersection[k]+fac*vecVT[j][k];
	      PRINTCOORD(c);
	      for (int k=0;k<NDIM;++k) c[k]+=fac*vecVN[j][k];
	      PRINTCOORD(c);
	      int index=0;if (dim==0) index++;if (j>0) index++;
	      if (dim==index) index++;
	      c[index]+=fac*vEdgeContribSign[index];
	      cs[j][index] = vEdgeContribSign[index];
	      PRINTCOORD(c);

	      for (int k=0;k<NDIM;++k) c[k]=intersection[k]+fac*vecVT[j][k];
	      PRINTCOORD(c);
	      for (int k=0;k<NDIM;++k) c[k]+=fac*vecSN[j][k];
	      PRINTCOORD(c);
	      for (int k=0;k<NDIM;++k) c[k]+=fac*facetNormal[k];
	      PRINTCOORD(c);
	    }
	  //fprintf(fl[0],"2 1\n 0 1 2\n[ADDITIONAL_DATA]\n");
	  fprintf(fl[0],"3 1\n 0 1 2 3\n[ADDITIONAL_DATA]\n");
	  for (int k=0;k<NDIM;++k)
	    {
	      if (k==0) fprintf(fl[0],"vecX\n0\n");
	      else if (k==1) fprintf(fl[0],"vecY\n0\n");
	      else fprintf(fl[0],"vecZ\n0\n");
	      fprintf(fl[0],"0 0 0 0\n");
	      fprintf(fl[0],"%e %e %e 0\n",fac*vecVT[0][k],fac*vecVN[0][k],fac*cs[0][k]);
	      fprintf(fl[0],"%e %e 0\n",fac*vecSN[0][k],fac*facetNormal[k]);
	      fprintf(fl[0],"%e %e %e 0\n",fac*vecVT[1][k],fac*vecVN[1][k],fac*cs[1][k]);
	      fprintf(fl[0],"%e %e 0\n",fac*vecSN[1][k],fac*facetNormal[k]);
	    }
	}
#endif	    
      return contrib;
      //return nProcessed;
    }

    // Correct the gradient contribution to deltaWeight if necessary, when an edge is 
    // crossing a periodic boundary
    template <class T, class WF>
    Float getPeriodizedDeltaWeight(const T  * coord, const WF &wf, 
				   const IncidentFace &face)
    {     
      Float result=face.deltaWeight; //return result;
      
      if ((AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)&&(WF::ORDER>0))
	{
	  // We need to be carefull here, as the gradient contribution is computed
	  // relative to the vertex with index 0 reference, so this vertex has to
	  // be on the same side as coord even when the simplex crosses a periodic
	  // boundary. If this is not the case, then we have to fix it !

	  const Coord* nextRefC = face.nextSimplex->getVertex(0)->getCoordsConstPtr();
	  const Coord* prevRefC = nextRefC;

	  // By convention, only prevSimplex may be NULL on a boundary !
	  if (face.prevSimplex!=NULL)
	    prevRefC=face.prevSimplex->getVertex(0)->getCoordsConstPtr();

	  T nextRef[3];
	  T prevRef[3];	  
	  std::copy_n(nextRefC,3,nextRef);
	  std::copy_n(prevRefC,3,prevRef);	 

	  T u0 = meshGeometry->checkCoordConsistency(nextRef[0],coord[0],0);
	  T u1 = meshGeometry->checkCoordConsistency(nextRef[1],coord[1],1);
	  T u2 = meshGeometry->checkCoordConsistency(nextRef[2],coord[2],2);
	  T v0 = meshGeometry->checkCoordConsistency(prevRef[0],coord[0],0);
	  T v1 = meshGeometry->checkCoordConsistency(prevRef[1],coord[1],1);
	  T v2 = meshGeometry->checkCoordConsistency(prevRef[2],coord[2],2);
	      
	  if ((u0!=nextRef[0])||(u1!=nextRef[1])||(u2!=nextRef[2])||
	      (v0!=prevRef[0])||(v1!=prevRef[1])||(v2!=prevRef[2]))
	    {
	      // deltaWeight needs to be updated as this is where we store the 
	      // contribution of the reference point
	      /*
	      Float nextWeight=wf.get(face.nextSimplex);	      		 
	      Float nextGrad[3];
	      wf.getGradient(face.nextSimplex,nextGrad);

	      // We compensate for the error on the ref point coordinates by previously
	      // substracting it from the 0th order term.
	      nextWeight -= 
		nextGrad[0]*(u0-nextRef[0])+ 
		nextGrad[1]*(u1-nextRef[1])+
		nextGrad[2]*(u2-nextRef[2]);
	      */
	      T coords[3]={u0,u1,u2};
	      Float nextWeight=wf.template compute<T,Float>(face.nextSimplex,coords);

	      Float prevWeight=0;
	      if (face.prevSimplex!=NULL)
		{
		  // not a boundary edge
		  /*
		  prevWeight=wf.get(face.prevSimplex);
		  Float prevGrad[3];	      
		  wf.getGradient(face.prevSimplex,prevGrad);
		  prevWeight -= 
		    prevGrad[0]*(v0-prevRef[0])+ 
		    prevGrad[1]*(v1-prevRef[1])+
		    prevGrad[2]*(v2-prevRef[2]);
		  */
		  T coords[3]={v0,v1,v2};
		  prevWeight=wf.template compute<T,Float>(face.prevSimplex,coords);
		}

	      if (face.flags & IncidentFace::flag_boundary) 
		prevWeight = 0; // boundary face
	      else if (face.flags & IncidentFace::flag_fold) 
		prevWeight = -prevWeight; // fold face

	      if (face.flags & IncidentFace::flag_invalidPrev) 
		prevWeight = 0; // prev neighbor is a ghost
	      if (face.flags & IncidentFace::flag_invalidNext) 
		nextWeight = 0; // next neighbor is a ghost

	      result = nextWeight-prevWeight;
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

      if (validSimplex(simplex))
	{
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
	      T u2 = meshGeometry->checkCoordConsistency(ref[2],coord[2],2);
		      
	      if ((u0!=ref[0])||(u1!=ref[1])||(u2!=ref[2]))
		{
		  /*	      	      
		  Float grad[3];
		  wf.getGradient(simplex,grad);
		  
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
		  result -= 
		    grad[0]*(u0-ref[0]) + 
		    grad[1]*(u1-ref[1]) +
		    grad[2]*(u2-ref[2]);
#pragma GCC diagnostic pop
		  */
		  T coords[3]={u0,u1,u2};
		  result=wf.template compute<T,Float>(simplex,coords);
		}
	    }
	  return result;
	}
      return 0;
    }

    template <class T, class T2, class WF,class OutputIterator>
    int addVoxelCornerContrib(T *coords, 
			      T2 vertexNeighborDir[Voxel::NNEI][NDIM],
			      Voxel *neighbors[Voxel::NNEI],
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
      static const Float dimFactorInv = 6.0;
     
      Float contrib0=dimFactorInv*
	getPeriodizedWeight(coords,wf,simplex);//,coordsAreConsistent);	
      Float grad[NDIM];
      if (WF::ORDER>0) wf.getGradient(simplex,grad);

      T newCoords[NDIM];

      for (unsigned long n=0;n<(1<<NDIM);++n) 
	{
	  Float contrib=0;
	  Float PTPN=1.0F;
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
	      // each vector contributes twice to order 1 part,
	      // with a factor of 1/4, 2/4 and 3/4.
	      // => 2*(1/4+2/4+3/4) = 3
	      contrib *= 3;
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
  
  private:

    template <class T, class T2, class WF>
    void getPeriodizedWeightAndGrad(const T  * coord, T2 &w, T2 *grad, 
				    const WF &wf, Simplex *simplex,
				    bool coordsAreConsistent=false) const
    {
      if (validSimplex(simplex))
	{
	  w = wf.template get<T2>(simplex);	  
	  if (WF::ORDER>0)
	    {
	      wf.getGradient(simplex,grad);
	      if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)//&&(!coordsAreConsistent))
		{
		  //const Coord* ref = simplex->getVertex(0)->getCoordsConstPtr();
		  T ref[NDIM];
		  for (int i=0;i<NDIM;++i)
		    ref[i]=simplex->getVertex(0)->getCoord(i);

		  T u0 = meshGeometry->checkCoordConsistency(ref[0],coord[0],0);
		  T u1 = meshGeometry->checkCoordConsistency(ref[1],coord[1],1);
		  T u2 = meshGeometry->checkCoordConsistency(ref[2],coord[2],2);
		  if ((u0!=ref[0])||(u1!=ref[1])||(u2!=ref[2]))
		    {
		      /*
		      w -= 
			grad[0]*(u0-ref[0]) + 
			grad[1]*(u1-ref[1]) +
			grad[2]*(u2-ref[2]);
		      */
		      //Coord coords[3]={u0,u1,u2};
		      ref[0]=u0;
		      ref[1]=u1;
		      ref[2]=u2;
		      w=wf.template compute<T,T2>(simplex,ref);
		    }
		}
	    }
	}
      else
	{
	  w=0;
	  if (WF::ORDER>0) std::fill_n(grad,NDIM,0);
	}
    }

    // NOTE: prevSimplex may be NULL but NOT nextSimplex !
    template <class WF>
    void getFace(IncidentEdge &edge, 
		 Simplex *nextSimplex,
		 Simplex *prevSimplex, 
		 Vertex *faceVertex, Vertex *oppVertex,
		 const WF &wf, IncidentFace &face)
    {      
      Float nextWeight;
      Float nextGrad[NDIM];

      face.flags=0;

      if (validSimplex(nextSimplex))
	{
	  nextWeight = wf.template get<Float>(nextSimplex);
	  if (WF::ORDER>0) wf.getGradient(nextSimplex,nextGrad);
	}
      else
	{
	  // The contribution of this face will be computed by another process
	  // so we do not count it here
	  face.flags|=IncidentFace::flag_invalidNext;
	  nextWeight = 0;
	  if (WF::ORDER>0) std::fill_n(nextGrad,NDIM,0);	
	}

      // Gram-Schmidt for vecN
      // FIXME: this may be a problem if a vertex falls on a segment ???
      /* Float version */   // uncomment me 
      /*   
      Float dot;
      meshGeometry->template getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
							 faceVertex->getCoordsConstPtr(),
							 face.vecN);      
      dot = meshGeometry->template dot_noCheck<Float,Float,NDIM>(edge.vecT,face.vecN);
      for (int i=0;i<NDIM;++i) face.vecN[i]-=dot*edge.vecT[i];
      dot = meshGeometry->template normalize_noCheck<Float,NDIM>(face.vecN);

      // And cross product for vecB (NOT gram-schmidt as it would break when a simplex
      // is almost or exactly flat)
      
      face.vecB[0]=edge.vecT[1]*face.vecN[2]-edge.vecT[2]*face.vecN[1];
      face.vecB[1]=edge.vecT[2]*face.vecN[0]-edge.vecT[0]*face.vecN[2];
      face.vecB[2]=edge.vecT[0]*face.vecN[1]-edge.vecT[1]*face.vecN[0];          

      // now orient face.vecB toward the remaining vertex in the simplex: the one that 
      // does not belong to the face (i.e. oppVertex)
      Float vecB[NDIM];
      meshGeometry->template getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
							 oppVertex->getCoordsConstPtr(),
							 vecB);
      dot=meshGeometry->template dot_noCheck<Float,Float,NDIM>(face.vecB,vecB);
      if (dot<0)
	{
	  for (int i=0;i<NDIM;++i) 
	    face.vecB[i]=-face.vecB[i];
	}
      */
      FaceVectors<NDIM,MESH,Float,HFloat>::
	computeNB(face,edge,faceVertex,oppVertex,meshGeometry);

      /* HFloat version */  // uncomment me
      /*
      HFloat dotH;
      HFloat vecHN[NDIM];
      HFloat vecHT[NDIM];

      for (int i=0;i<NDIM;++i)
	vecHT[i]=hlp::numericStaticCast<HFloat>(edge.vecT[i]);

      // Gram-Schmidt for vecN
      // FIXME: this may be a problem if a vertex falls on a segment ???
      meshGeometry->template getVector<Coord,HFloat,NDIM>(edge.vertex->getCoordsConstPtr(),
							  faceVertex->getCoordsConstPtr(),
							  vecHN);  
      dotH = meshGeometry->template dot_noCheck<HFloat,HFloat,NDIM>(vecHT,vecHN);
      for (int i=0;i<NDIM;++i) vecHN[i]-=dotH*vecHT[i];
      dotH = meshGeometry->template normalize_noCheck<HFloat,NDIM>(vecHN);

      // And cross product for vecB (NOT gram-schmidt as it would break when a simplex
      // is almost or exactly flat)
      HFloat vecHB[NDIM];
      vecHB[0]=vecHT[1]*vecHN[2]-vecHT[2]*vecHN[1];
      vecHB[1]=vecHT[2]*vecHN[0]-vecHT[0]*vecHN[2];
      vecHB[2]=vecHT[0]*vecHN[1]-vecHT[1]*vecHN[0];
      
      // now orient face.vecB toward the remaining vertex in the simplex: the one that 
      // does not belong to the face (i.e. oppVertex)
      HFloat vecHB_dir[NDIM];
      meshGeometry->template getVector<Coord,HFloat,NDIM>(edge.vertex->getCoordsConstPtr(),
							  oppVertex->getCoordsConstPtr(),
							  vecHB_dir);
      dotH=meshGeometry->template dot_noCheck<HFloat,HFloat,NDIM>(vecHB,vecHB_dir); 
      */
      /* TEST */
      /*
      for (int i=0;i<NDIM;++i)
	{
	  HFloat hb=(dotH<0)?-vecHB[i]:vecHB[i];
	  if (hb != 0)
	    {
	      Float err=(hb-face.vecB[i]) / hb;
	      if (fabs(err) > 1.E-9) 
		std::cout << "vecB differ ("<<i<<","<<fabs(err)<<"): (" 
			  << vecHB[0] <<","<< vecHB[1] <<","<< vecHB[2]<<")("
			  << face.vecB[0] <<","<< face.vecB[1] <<","<< face.vecB[2]<<").\n";
	    }

	  if (vecHN[i] != 0)
	    {
	      Float err=(vecHN[i]-face.vecN[i]) / vecHN[i];
	      if (fabs(err) > 1.E-9) 
		std::cout << "vecN differ ("<<i<<","<<fabs(err)<<"): (" 
		  << vecHN[0] <<","<< vecHN[1] <<","<< vecHN[2]<<")("
		  << face.vecN[0] <<","<< face.vecN[1] <<","<< face.vecN[2]<<").\n";
	    }
	}
      */
      /* */

       // uncomment me
      /*
      if (dotH<0)
	{
	  for (int i=0;i<NDIM;++i)
	    {
	      face.vecB[i]=-hlp::numericStaticCast<Float>(vecHB[i]);
	      face.vecN[i]=hlp::numericStaticCast<Float>(vecHN[i]);
	    }
	}
      else
	{
	  for (int i=0;i<NDIM;++i)
	    {
	      face.vecB[i]=hlp::numericStaticCast<Float>(vecHB[i]);
	      face.vecN[i]=hlp::numericStaticCast<Float>(vecHN[i]);
	    }
	}
      */
      
      
      // if (dotH<0)
      // 	{
      // 	  for (int i=0;i<NDIM;++i) 
      // 	    face.vecB[i]=-face.vecB[i];
      // 	}
      
      /*
      // Then vecB. We don't actually need to do that for vecB but this way ensures a 
      // better othogonality with vecN and vecT ...
      dot = meshGeometry->template dot_noCheck<Coord,NDIM>(edge.vecT,face.vecB);
      for (int i=0;i<NDIM;++i) face.vecB[i]-=dot*edge.vecT[i];      
      dot = meshGeometry->template dot_noCheck<Coord,NDIM>(face.vecN,face.vecB);
      for (int i=0;i<NDIM;++i) face.vecB[i]-=dot*face.vecN[i];      
      dot = meshGeometry->template normalize_noCheck<Coord,NDIM>(face.vecB);
      */

      /*
      if (dot==0)
	{
	  // The binormal is degenerate => the simplex is flat !
	  // We use the facet normal instead, which will work for a first order
	  // degeneracy (i.e. tetrahedron reduces to a triangle, but not a segment/point)
	  int facetId;
	  for (int j=0;j<Simplex::NVERT;++j)
	    {
	      Vertex *v=nextSimplex->getVertex(j);
	      if ((v!=faceVertex)&&
		  (v!=edge.vertex)&&
		  (v!=edge.otherVertex))
		facetId=j;
	    }

	  // FIXME: don't we have to check the normal orientation here ?
	  nextSimplex->getFacetHandle(facetId)->
	    computeProjectedNormal(face.vecB,meshGeometry);
	  dot = meshGeometry->template normalize_noCheck<Coord,NDIM>(face.vecB);

	  face.vecN[0]=edge.vecT[1]*face.vecB[2]-edge.vecT[2]*face.vecB[1];
	  face.vecN[1]=edge.vecT[2]*face.vecB[0]-edge.vecT[0]*face.vecB[2];
	  face.vecN[2]=edge.vecT[0]*face.vecB[1]-edge.vecT[1]*face.vecB[0];

	  if (dot==0)
	    {	      
	      PRINT_SRC_INFO(LOG_WARNING);
	      glb::console->print<LOG_WARNING>("I encountered a highly degenerate tetrahedron (reduces to a segment or point). Cross your fingers ;)\n");
	      edge.vecT[0]=1;edge.vecT[1]=0;edge.vecT[2]=0;
	      face.vecN[0]=0;face.vecN[1]=1;face.vecN[2]=0;
	      face.vecB[0]=0;face.vecB[1]=0;face.vecB[2]=1;
	    }	 
	}
      */
      //if (validSimplex(prevSimplex))
      if (prevSimplex!=NULL)
	{
	  //face.prevWeight = prevSimplex->cache.d;
	  //face.prevWeight = validSimplex(prevSimplex)?prevSimplex->cache.d:0;
	  if (validSimplex(prevSimplex))
	    {
	      face.deltaWeight = wf.template get<Float>(prevSimplex);
	      if (WF::ORDER>0) wf.getGradient(prevSimplex,face.deltaGrad);
	    }
	  else
	    {
	      face.flags|=IncidentFace::flag_invalidPrev;
	      face.deltaWeight = 0;
	      if (WF::ORDER>0) std::fill_n(face.deltaGrad,NDIM,0);	
	    }

	  // We have to check that the facet between prevSimplex and nextSimplex
	  // is not a fold (caustic)...
	  // oppVertex2 will store the equivalent of oppVertex but for the prevSimplex
	  Vertex *oppVertex2;
	  for (int j=0;j<Simplex::NVERT;++j)
	    {
	      oppVertex2 = prevSimplex->getVertex(j);
	      if ((oppVertex2!=faceVertex)&&
		  (oppVertex2!=edge.vertex)&&
		  (oppVertex2!=edge.otherVertex))
		break;	
	    }

	  Float tmpVec[NDIM];
	  meshGeometry->template
	    getVector<Coord,Float,NDIM>(edge.vertex->getCoordsConstPtr(),
					oppVertex2->getCoordsConstPtr(),
					tmpVec);

	  // If it is a fold, then its contribution should be negated because
	  // its binormal will point toward the same direction as that of nextSimplex
	  Float dot=meshGeometry->template dot_noCheck<Float,Float,NDIM>(face.vecB,tmpVec);
	  //if (dot>0) {face.prevWeight = -face.prevWeight;}
	  if (dot>0) 
	    {
	      // fold face
	      face.deltaWeight = -face.deltaWeight;
	      if (WF::ORDER>0)
		{		  
		  face.deltaGrad[0] = -face.deltaGrad[0];
		  face.deltaGrad[1] = -face.deltaGrad[1];
		  face.deltaGrad[2] = -face.deltaGrad[2];
		  face.flags|=IncidentFace::flag_fold;
		}
	    }
	}
      else 
	{
	  // boundary face
	  //face.prevWeight = 0; 
	  face.deltaWeight = 0;
	  if (WF::ORDER>0)
	    {
	      std::fill_n(face.deltaGrad,NDIM,0);	      
	      face.flags|=IncidentFace::flag_boundary; // boundary face
	    }
	}

      // FIXME we could optimize here by skipping when both are 0 ...      
      //face.nextWeight = validSimplex(nextSimplex)?nextSimplex->cache.d:0;
      
      face.deltaWeight = nextWeight - face.deltaWeight;
      /*
      if (face.deltaWeight==0) face.deltaWeight = nextWeight;
      else if (nextWeight==0) face.deltaWeight = - face.deltaWeight;      
      else face.deltaWeight = wf.getDeltaWeight(nexSimplex,prevSimplex);
      */

      if (WF::ORDER>0)
	{		  
	  face.deltaGrad[0] = nextGrad[0] - face.deltaGrad[0];
	  face.deltaGrad[1] = nextGrad[1] - face.deltaGrad[1];
	  face.deltaGrad[2] = nextGrad[2] - face.deltaGrad[2];
	}
    
      // face.prevWeight = validSimplex(prevSimplex)?prevSimplex->cache.d:0;
      // face.nextWeight = validSimplex(nextSimplex)?nextSimplex->cache.d:0;

      face.prevSimplex = prevSimplex;
      face.nextSimplex = nextSimplex;
    }

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

    template <class T>
    bool findVertex(double x1,double y1,double z1,const T *vc, double tol=2.E-5) const
    {
      if ((fabs(vc[0]-x1) <= tol)&&
	  (fabs(vc[1]-y1) <= tol)&&
	  (fabs(vc[2]-z1) <= tol))
	{
	  return true;
	}
      return false;
    }
    
#ifdef DEBUGDUMP
    int nFiles;
    FILE *fl[10];
    void openFile(const char *name)
    {
      static int count[10]={0,0,0,0,0,0,0,0,0};
      char fname[255];
      sprintf(fname,"%s_%6.6d.dat",name,count[nFiles]++);
      fl[nFiles++]=fopen(fname,"w");
    }

    void closeFile()
    {
      for (int i=0;i<nFiles;++i) fclose(fl[i]);
      nFiles=0;
    }  
#endif

  };

}

#ifdef DEBUGDUMP
#undef PRINTCOORD
#undef PRINTVERTEX
#undef DEBUGDUMP
#endif

#ifdef DEBUG_CHECK
#undef DEBUG_CHECK_RANK
#undef DEBUG_CHECK
#endif

#include "../../internal/namespace.footer"
#endif
