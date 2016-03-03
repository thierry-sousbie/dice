#ifndef __FACET_HXX__
#define __FACET_HXX__

#include <algorithm>

#include "../mesh/internal/getSegmentHandles_internal.hxx"

#include "../mesh/vertex.hxx"
#include "../mesh/segment.hxx"
#include "../mesh/simplex.hxx"

#include "../mesh/internal/facetImplementation.hxx"

#include "../geometry/predicates/pointInSimplex.hxx"

/**
 * @file 
 * @brief Defines a facet class to be used with MeshT (used to store any type of cell)
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

/**
 * \class FacetT
 * \brief Defines a facet class to be used with MeshT. Facets are designed to be able to
 * store any dimension cells (i.e. simplices or any of their k-faces)
 */

template <class T>
class FacetT
{
private:
  typedef internal::FacetImplementationT< T::NDIM > Implementation;

public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class SimplexT<T>;
  friend class SegmentT<T>;
  friend class VertexT<T>;
  
  static const long NDIM = T::NDIM;
  static const long NDIM_W = T::NDIM_W;
  static const long NVERT = NDIM;
  static const long NSEG = (NDIM*(NDIM-1))/2;

  typedef FacetT<T> MyType;
  typedef T Traits;  
  
  typedef SimplexT<T> Simplex;
  typedef SegmentT<T> Segment;
  typedef FacetT<T> Facet;
  typedef VertexT<T> Vertex;

  typedef typename Vertex::Coord Coord;

  typedef typename Vertex::GlobalIdentity GlobalIdentity;
  typedef typename Vertex::LocalIndex LocalIndex;

  typedef typename Simplex::FacetHandle        Handle;
  typedef typename Simplex::SegmentHandle      SegmentHandle;
  typedef typename Simplex::ConstFacetHandle   ConstHandle;
  typedef typename Simplex::ConstSegmentHandle ConstSegmentHandle;
  
  FacetT(Simplex *s, int vIndex)
  {
    simplex=s;
    vertexIndex=vIndex;
  }

  FacetT()
  {
    simplex = NULL;
    vertexIndex=-1;
  }

  bool isSet() const
  {
    return ((simplex != NULL)&&(vertexIndex>=0));
  }

  Vertex *getOppositeVertex() const
  {
    return simplex->getVertex(vertexIndex);
  }

  char getOppositeVertexIndex() const
  {
    return vertexIndex;
  }
  
  Vertex *getVertex(int which) const
  {
    return Implementation::template getVertex<>(simplex,vertexIndex,which);
    //return (which<vertexIndex)?simplex->getVertex(which):simplex->getVertex(which+1);
  }
  
  template <class OutputIterator>
  void getVertices(OutputIterator out) const
  {
    Implementation::template getVertices<>(simplex,vertexIndex,out);  
  }

  template <typename OutputIterator>
  void getVerticesCoordsPtr(OutputIterator coords) const
  {
    Vertex *tmp[NVERT];
    getVertices(tmp);
    for (long i=0;i<NVERT;i++)
      {
	(*coords)=tmp[i]->getCoordsPtr();
	++coords;
      }
  }

  template <typename OT>
  void getVerticesCoordsPtr(OT * /*__restrict*/ coords) const
  {
    Vertex *tmp[NVERT];
    getVertices(tmp);
    for (long i=0;i<NVERT;i++)
      coords[i]=tmp[i]->getCoordsPtr();   
  }

  template <typename OutputIterator>
  void getVerticesCoordsConstPtr(OutputIterator coords) const
  {
    Vertex *tmp[NVERT];
    getVertices(tmp);
    for (long i=0;i<NVERT;i++)
      {
	(*coords)=tmp[i]->getCoordsConstPtr();
	++coords;
      }
  }

  template <typename OT>
  void getVerticesCoordsConstPtr(OT* /*__restrict*/ coords) const
  {
    Vertex *tmp[NVERT];
    getVertices(tmp);
    for (long i=0;i<NVERT;i++)
      coords[i]=tmp[i]->getCoordsConstPtr();     
  }

  template <class OutputIterator>
  void getVerticesLocalIndex(OutputIterator out)  const
  {
    Vertex *tmp[NVERT];
    getVertices(tmp);
    for (int i=0;i<NVERT;++i)
      {
	*out=tmp[i]->getLocalIndex();
	++out;
      }
  }
  

  template <class OutputIterator, class G, class CT=Coord>
  void computeProjectedNormal(OutputIterator out, const G *geometry) const
  {
    CT vCoord[NVERT][NDIM];

    for (int n=0;n<NVERT;++n)
      {
	const Coord *c=getVertex(n)->getCoordsConstPtr();
	for (int d=0;d<NDIM;++d)
	  vCoord[n][d]=c[d];
      }    
    
    geometry->template checkCoordsConsistency<CT,NVERT,NDIM>(vCoord,vCoord[0]);

    Implementation::template computeProjectedNormal<CT,OutputIterator>(vCoord,out);
  }

  template <class OutputIterator, class CT=Coord>
  void computeProjectedNormal(OutputIterator out) const
  {
    CT vCoord[NVERT][NDIM];

    for (int n=0;n<NVERT;++n)
      {
	const Coord *c=getVertex(n)->getCoordsConstPtr();
	for (int d=0;d<NDIM;++d)
	  vCoord[n][d]=c[d];
      }
    
    Implementation::template computeProjectedNormal<CT,OutputIterator>(vCoord,out);
  }

  SegmentHandle getSegmentHandle(int which)  const
  {
    Vertex *vertices[NVERT];
    getVertices(vertices);

    return internal::GetSegmentHandle<NDIM-1,SegmentHandle,Segment>::
      get(which,simplex,vertices);
  }

  template <class OutputIterator>
  void getSegmentHandles(OutputIterator out)  const
  {
    Vertex *vertices[NVERT];
    getVertices(vertices);

    for (int i=0;i<NSEG;++i)
      {
	(*out) = internal::GetSegmentHandle<NDIM-1,SegmentHandle,Segment>::
	  get(i,simplex,vertices);
	++out;
      }
  }

  Simplex *getSimplex() const 
  {
    return simplex;
  }
  
  Simplex *getOppositeSimplex() const 
  {
    return simplex->getNeighbor(vertexIndex);
  }
  
  Simplex *swapSimplex()
  {
    Simplex *other = getOppositeSimplex();

    if (other==NULL) return NULL;

    Vertex *ov[Simplex::NVERT];
    Vertex *v[NVERT];

    other->getVertices(ov);
    getVertices(v);

    int i=0;
    for (;i<Simplex::NVERT;++i)
      {	
	bool found=false;
	for (int j=0;j<NVERT;++j)
	  if (v[j]==ov[i]) found=true;
	if (!found) break;
      }

    vertexIndex = i;
    simplex = other;
    
    return simplex;
  }
  
  Handle getHandle()
  {
    return Handle(this);
  }

  ConstHandle getConstHandle() const
  {
    return ConstHandle(this);
  }

  template <class L>
  void print() const
  {
    glb::console->print<L>("Facet @%d from simplex :\n",(int)vertexIndex); 
    simplex->template print<L>();
    for (int i=0;i<NVERT;++i)
      getVertex(i)->template print<L>();  
  }
  /*
  template <class G, class TZ>
  int intersectRay(const Coord * coords, const G * geometry,
		   int dim, TZ &intersectionCoord) const
  {
    return predicate::PointInSimplexT<NDIM>::
      template intersectRayFacet<Simplex,G,TZ>
      (simplex,coords,geometry,dim,vertexIndex,intersectionCoord); 
  }
  */
  template <class G, class TZ, class TC=TZ, int filter = predicate::filterType::Adaptive>
  int intersectSegment(const Coord * coords, int dim, Coord otherCoord,
		       const G * geometry, TZ &intersectionCoord, 
		       int coordsAreConsistent=-1) const
  {
    return predicate::PointInSimplexT<NDIM,filter>::
      template intersectSegmentFacet<MyType,G,TZ,TC>
      (*this,coords,dim,otherCoord,geometry,intersectionCoord,coordsAreConsistent); 
  }
  
protected:
  Simplex *simplex;
  char vertexIndex;
};

/** \}*/
#include "../internal/namespace.footer"
#endif
