#ifndef __SIMPLEX_FROM_VERTICES_HXX__
#define __SIMPLEX_FROM_VERTICES_HXX__

#include <string.h> //for memset

#include <vector>
#include <algorithm>

#include "../mesh/internal/getSegmentHandles_internal.hxx"

#include "../mesh/vertex.hxx"
#include "../mesh/segment.hxx"
#include "../mesh/facet.hxx"
#include "../tools/helpers/helpers.hxx"
#include "../tools/types/handle.hxx"

/**
 * @file 
 * @brief Defines an implementation for SimplexT that store each vertices of the simplex.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

template <class T> class LocalMeshT;
template <class T> class MeshT;

/**
 * \class SimplexFromVerticesT
 * \brief Defines an implementation for SimplexT that store each vertices of the simplex.
 */

template <class T, template <class> class S>
class SimplexFromVerticesT : public T::SimplexData
{
public:
  static const long NDIM     = T::NDIM;
  static const long NDIM_W   = T::NDIM_W;

  static const long NVERT    = NDIM+1;
  static const long NSEG     = (NDIM*(NDIM+1))/2;
  static const long NFACET   = NDIM+1;
  static const long NNEI     = NVERT;
  static const long NOPPSEG  = NSEG-NDIM;
  static const long NOPPVERT = NVERT-1;
  static const long NOPPEL   = NOPPVERT;  

  static const bool EXPLICIT_SEGMENTS = false;  

  typedef typename T::SimplexFlag    Flag;
  typedef typename T::SimplexData    Data;
  typedef typename T::Coord          Coord;
  typedef typename T::LocalIndex     LocalIndex;
  typedef typename T::GlobalIndex    GlobalIndex;
  typedef typename T::GlobalIdentity GlobalIdentity;

  static const LocalIndex  LOCAL_INDEX_INVALID;//  = T::LOCAL_INDEX_INVALID;
  static const LocalIndex  LOCAL_INDEX_MAX;//      = T::LOCAL_INDEX_MAX;
  static const GlobalIndex GLOBAL_INDEX_INVALID;// = T::GLOBAL_INDEX_INVALID;
  static const GlobalIndex GLOBAL_INDEX_MAX;//     = T::GLOBAL_INDEX_MAX;
 
  typedef SimplexFromVerticesT<T,S> MyType;  
 
  typedef S<T>        Simplex; 
  typedef FacetT<T>   Facet;   /**< see FacetT */
  typedef SegmentT<T> Segment; /**< see SegmentT */
  typedef VertexT<T>  Vertex;  /**< see VertexT */
  typedef Vertex      Element;

  typedef HandleT<Segment>      SegmentHandle;
  typedef HandleT<Facet>        FacetHandle;
  typedef ConstHandleT<Segment> ConstSegmentHandle;
  typedef ConstHandleT<Facet>   ConstFacetHandle;

  SimplexFromVerticesT():
    localIndex(0),
    flags(SIMPLEX_FLAG_NOTSET|SIMPLEX_FLAG_UNSAFE)    
  {
    Simplex* nullSimplex = NULL;
    Vertex* nullVertex = NULL;
    //Segment* nullSegment = NULL;
    std::fill(neighbors,neighbors+NNEI,nullSimplex);
    std::fill(vertices,vertices+NVERT,nullVertex);
    //std::fill(segments,segments+NSEG,nullSegment);
  }

  ~SimplexFromVerticesT() {}

  // Used when a simplex is refined
  void shrink(const SegmentHandle &seg, Vertex *newVertex, Segment **newSegment, 
	      bool direct, Simplex *partner)
  {
    Vertex *v0=seg->getVertex(0);
    Vertex *v1=seg->getVertex(1);
    int i0=-1;
    int i1=-1;
    
    // ensure that whatever the segment, the ordering of the split simplices 
    // is always the same
    if ((*v0)>(*v1)) std::swap(v0,v1); 

    for (int i=0;i<NVERT;i++)
      {
	if (v0==vertices[i]) i0=i;
	if (v1==vertices[i]) i1=i;
	//if (neighbors[i]==NULL) printf("NULL\n");
      }

    if ((i0<0)||(i1<0))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Invalid simplex or segment handle.\n");	    
	exit(-1);
      }
    
    // If 'direct=true', the segment's vertex with lowest index in vertices is preserved 
    // and the other is moved, the opposite happens when 'direct=false'.
    // This ensures consistency over the split simplices vertices order, and also
    // preserves simplices orientation !
    int dir=((direct)?1:-1)*((i0<i1)?1:-1);
    // (*v) is where we will put the new vertex
    Vertex **v = (dir>0)?(&vertices[i1]):(&vertices[i0]);
    // (*s) is where the other part of the split simplex will be neighbor
    Simplex **s = (dir<0)?(&neighbors[i1]):(&neighbors[i0]);
    
    *v = newVertex;
    *s = partner; 
     
    setSafeF(false);
  }
  
  int findSegmentIndex(const SegmentHandle &seg) const
  {
    return findSegmentIndex(seg->getVertex(0),seg->getVertex(1));
  }

  int findSegmentIndex(const ConstSegmentHandle &seg) const
  {
    return findSegmentIndex(seg->getVertex(0),seg->getVertex(1));
  }

  int findSegmentIndex(const Vertex *v1,const Vertex *v2) const
  {
    int result=-1;
    int i1=getVertexIndex(v1);
    int i2=getVertexIndex(v2);

    if ((i1>=0)&&(i2>=0))
      result=findSegmentIndex(i1,i2);
    
    return result;
  }

  int findSegmentIndex(int vid1, int vid2) const
  {
    if (vid1>vid2) std::swap(vid1,vid2);

    if (vid1==vid2) 
      return -1;
    else 
      return NDIM*vid1 - (vid1*(vid1-1))/2 + vid2-vid1-1;    
  }

  // vid1 must me less than vid2
  template <int VID1, int VID2>
  static int findSegmentIndex()
  {
    static const int vid1= hlp::MinT<VID1,VID2>::value;
    static const int vid2= hlp::MaxT<VID1,VID2>::value;
    return NDIM*vid1 - (vid1*(vid1-1))/2 + vid2-vid1-1;    
  }
  
  /*
  int findSegmentIndex(int vid1, int vid2) const
  {    
    if (vid1==vid2) return -1;
    if (vid1>vid2) std::swap(vid1,vid2);
    
    int j=0;
    int index=(vid2-vid1-1);

    while (j!=vid1)
      {
	index+=(NDIM-j);
	j++;
      };

    return index;
  }
  */

  SegmentHandle getSegmentHandle(int i)
  {  
    return internal::GetSegmentHandle<NDIM,SegmentHandle,Segment>::
      get(i,static_cast<Simplex*>(this),vertices);    
  }


  ConstSegmentHandle getConstSegmentHandle(int i) const
  {
    return internal::GetSegmentHandle<NDIM,ConstSegmentHandle,Segment>::
      get(i,static_cast<const Simplex*>(this),vertices);
  }

  static std::pair<int,int> getSegmentVerticesIndex(int i)
  {
    return internal::GetSegmentHandle<NDIM,SegmentHandle,Segment>::
      getIndices(i);
  }
 
  FacetHandle getFacetHandle(int i) 
  {
    return FacetHandle(Facet(static_cast<Simplex *>(this),i));
  }

  ConstFacetHandle getConstFacetHandle(int i) const
  {
    return ConstFacetHandle(Facet(static_cast<Simplex*>(const_cast<MyType *>(this)),i));
  }
  
  template <typename OutputIterator>
  void getVerticesCoordsPtr(OutputIterator coords) const
  {
    for (long i=0;i<NVERT;i++)
      {
	(*coords)=vertices[i]->getCoordsPtr();
	++coords;
      }
  }

  template <typename OT>
  void getVerticesCoordsPtr(OT * /*__restrict*/ coords) const
  {
    for (long i=0;i<NVERT;i++)
      coords[i]=vertices[i]->getCoordsPtr();   
  }

  template <typename OutputIterator>
  void getVerticesCoordsConstPtr(OutputIterator coords) const
  {
    for (long i=0;i<NVERT;i++)
      {
	(*coords)=vertices[i]->getCoordsConstPtr();
	++coords;
      }
  }

  template <typename OT>
  void getVerticesCoordsConstPtr(OT* /*__restrict*/ coords) const
  {
    for (long i=0;i<NVERT;i++)
      coords[i]=vertices[i]->getCoordsConstPtr();     
  }
  
  template <typename OutputIterator>
  void getVertices(OutputIterator vert) const
  {
    for (int i=0;i<NVERT;++i)
      {
	(*vert)=vertices[i];
	++vert;
      }
  }
  
  //template <typename OutputIterator>
  void getVertices_restrict(Vertex ** __restrict vert) const
  {
    Vertex *const* __restrict cur=vertices;
    std::copy(cur,cur+NVERT,vert);
    /*
    for (int i=0;i<NVERT;++i)
      vert[i]=cur[i];
    */
  }

  template <typename OutputIterator>
  void getVerticesLocalIndex(OutputIterator vert) const
  {
    for (long i=0;i<NVERT;i++)
      {
	(*vert)=vertices[i]->getLocalIndex();
	++vert;
      }    
  }

  int getNeighborIndex(const MyType *nei)
  {
    for (int i=0;i<NNEI;++i)
      if (neighbors[i]==nei) return i;
    return -1;
  }
  
  Vertex *getVertex(int id) const
  {
    return vertices[id];
  }

  Vertex *getVertexByGlobalIdentity(typename Vertex::GlobalIdentity id) const
  {
    for (int i=0;i<NVERT;i++)
      if (vertices[i]->getGlobalIdentity() == id) return vertices[i];
    return NULL;
  }

  Vertex *getVertexByLocalIndex(typename Vertex::LocalIndex id) const
  {
    for (int i=0;i<NVERT;i++)
      if (vertices[i]->getLocalIndex() == id) return vertices[i];
    return NULL;
  }
  
  int getVertexIndex(const Vertex *v) const
  {   
    for (long i=0;i<NVERT;i++)    
      if (vertices[i]==v) return i;
      
    return -1;
  }
  
  // get all the vertices except 'ref'
  template <typename OutputIterator>
  void getOppositeVertices(OutputIterator vert, int refVertex) const
  {
    for (long i=0;i<NVERT;i++)
      {
	if (i!=refVertex)
	  {
	    (*vert)=vertices[i];
	    ++vert;
	  }
      }   
  }  

  template <typename OutputIterator>
  void getOppositeElements(OutputIterator vert, int refVertex) const
  {
    getOppositeVertices(vert,refVertex);
  }

  const Data &getData() const
  {
    return *static_cast<const Data*>(this);
  }

  Data &getData()
  {
    return *static_cast<Data*>(this);
  }

  void setData(const Data &data)
  {
    Data &d=getData();
    d=data;
  }

  Flag getFlags() const {return flags;}

  bool isSafe() const
  {
    return !(flags & SIMPLEX_FLAG_UNSAFE);
  }

  bool isSet() const
  {
    return !(flags & SIMPLEX_FLAG_NOTSET);
  }
  
  bool isShared() const
  {
    return flags & SIMPLEX_FLAG_SHARED;
  }
  
  bool isBoundary() const
  {
    return flags & SIMPLEX_FLAG_BOUNDARY;
  }

  bool isShadow() const
  {
    return flags & SIMPLEX_FLAG_SHADOW;
  }

  bool isGhost() const
  {
    return flags & SIMPLEX_FLAG_GHOST;
  }

  bool isShadowOrGhost() const
  {
    return flags & (SIMPLEX_FLAG_SHADOW|SIMPLEX_FLAG_GHOST);
  }

  bool isTagged() const
  {
    return (flags&SIMPLEX_FLAG_TAG);
  }

  bool isLocal() const
  {
    return !isShadowOrGhost();
  }

  template <class L>
  void print(const std::string &s) const
  {
    if (NDIM==1)
      {
//#if (DEFINED_MESH_DIM==1) //this is to prevent spurious warnings from icc
	glb::console->print<L>("%sSimplex %ld(%ld) V=[%ld,%ld], F=%d, Nei=[%ld(%d),%ld(%d)]\n",
			  s.c_str(),(long)static_cast<const Simplex*>(this),
			  (long)localIndex,
			  (long)vertices[0],(long)vertices[1],
			  flags,
			  (long)neighbors[0],
			  (neighbors[0]==NULL)?-1:(int)neighbors[0]->flags,
			  (long)neighbors[1],
			  (neighbors[1]==NULL)?-1:(int)neighbors[1]->flags);
//#endif // (DEFINED_MESH_==1)
      }
    else if (NDIM==2)
      {
//#if (DEFINED_MESH_DIM==2) //this is to prevent spurious warnings from icc
	glb::console->print<L>("%sSimplex %ld(%ld) V=[%ld(%d),%ld(%d),%ld(%d)]=[(%ld,%ld);(%ld,%ld);(%ld,%ld)], F=%d, Nei=[%ld(%d),%ld(%d),%ld(%d)]=[(%ld,%ld);(%ld,%ld);(%ld,%ld)]\n",
			  s.c_str(),(long)static_cast<const Simplex*>(this),
			  (long)localIndex,
			  (long)vertices[0],
			  (vertices[0]==NULL)?-1:(int)vertices[0]->getFlags(),
			  (long)vertices[1],
			  (vertices[1]==NULL)?-1:(int)vertices[1]->getFlags(),
			  (long)vertices[2],
			  (vertices[2]==NULL)?-1:(int)vertices[2]->getFlags(),
			  (vertices[0]==NULL)?-1:(long)vertices[0]->getGlobalIdentity().rank(),
			  (vertices[0]==NULL)?-1:(long)vertices[0]->getGlobalIdentity().id(),
			  (vertices[1]==NULL)?-1:(long)vertices[1]->getGlobalIdentity().rank(),
			  (vertices[1]==NULL)?-1:(long)vertices[1]->getGlobalIdentity().id(),
			  (vertices[2]==NULL)?-1:(long)vertices[2]->getGlobalIdentity().rank(),
			  (vertices[2]==NULL)?-1:(long)vertices[2]->getGlobalIdentity().id(),
			  flags,
			  (long)neighbors[0],
			  (neighbors[0]==NULL)?-1:(int)neighbors[0]->flags,
			  (long)neighbors[1],
			  (neighbors[1]==NULL)?-1:(int)neighbors[1]->flags,
			  (long)neighbors[2],
			  (neighbors[2]==NULL)?-1:(int)neighbors[2]->flags,
			  (neighbors[0]==NULL)?-1:(long)neighbors[0]->getGlobalIdentity(GlobalIdentity::MAX_RANK).rank(),
			  (neighbors[0]==NULL)?-1:(long)neighbors[0]->getGlobalIdentity(GlobalIdentity::MAX_RANK).id(),
			  (neighbors[1]==NULL)?-1:(long)neighbors[1]->getGlobalIdentity(GlobalIdentity::MAX_RANK).rank(),
			  (neighbors[1]==NULL)?-1:(long)neighbors[1]->getGlobalIdentity(GlobalIdentity::MAX_RANK).id(),
			  (neighbors[2]==NULL)?-1:(long)neighbors[2]->getGlobalIdentity(GlobalIdentity::MAX_RANK).rank(),
			  (neighbors[2]==NULL)?-1:(long)neighbors[2]->getGlobalIdentity(GlobalIdentity::MAX_RANK).id());
//#endif // (DEFINED_MESH_DIM==2)
      }
    else if (NDIM==3)
      {
//#if (DEFINED_MESH_DIM==3) //this is to prevent spurious warnings from icc
	glb::console->
	  print<L>("%sSimplex %ld::%d V=[%ld,%ld,%ld,%ld]=[(%ld,%ld);(%ld,%ld);(%ld,%ld);(%ld,%ld)], F=%d, Nei=[%ld(%d),%ld(%d),%ld(%d),%ld(%d)]\n",
		   s.c_str(),(long)static_cast<const Simplex*>(this),
		   (long)localIndex,			  
		   (long)vertices[0],
		   (long)vertices[1],
		   (long)vertices[2],
		   (long)vertices[3],
		   (vertices[0]==NULL)?-1:(long)vertices[0]->getGlobalIdentity().rank(),
		   (vertices[0]==NULL)?-1:(long)vertices[0]->getGlobalIdentity().id(),
		   (vertices[1]==NULL)?-1:(long)vertices[1]->getGlobalIdentity().rank(),
		   (vertices[1]==NULL)?-1:(long)vertices[1]->getGlobalIdentity().id(),
		   (vertices[2]==NULL)?-1:(long)vertices[2]->getGlobalIdentity().rank(),
		   (vertices[2]==NULL)?-1:(long)vertices[2]->getGlobalIdentity().id(),
		   (vertices[3]==NULL)?-1:(long)vertices[3]->getGlobalIdentity().rank(),
		   (vertices[3]==NULL)?-1:(long)vertices[3]->getGlobalIdentity().id(),
		   flags,
		   (long)neighbors[0],
		   (neighbors[0]==NULL)?-1:(int)neighbors[0]->flags,
		   (long)neighbors[1],
		   (neighbors[1]==NULL)?-1:(int)neighbors[1]->flags,
		   (long)neighbors[2],
		   (neighbors[2]==NULL)?-1:(int)neighbors[2]->flags,
		   (long)neighbors[3],
		   (neighbors[3]==NULL)?-1:(int)neighbors[3]->flags);
//#endif // (DEFINED_MESH_DIM==3)
      }
  }

  template <class W>
  static void selfSerialize(const MyType *me, W *writer)
  {
    Data::selfSerialize(static_cast<const Data*>(me), writer);
    writer->write(me->vertices,NVERT);
    writer->write(me->neighbors,NNEI);
    writer->write(&me->generation);
    writer->write(&me->localIndex);
    writer->write(&me->flags);        
  }

  template <class R>
  static void selfUnSerialize(MyType *me, R *reader)
  {
    Data::selfUnSerialize(static_cast<Data*>(me), reader);
    reader->read(me->vertices,NVERT);
    reader->read(me->neighbors,NNEI);
    reader->read(&me->generation);
    reader->read(&me->localIndex);
    reader->read(&me->flags);    
  }

  union Cache
  {
    //std::pair<float,int> pfi;
    struct {float f; int i;} pfi;
    char c[8];    
    int i[2];
    unsigned int ui[2];
    float f[2];
    double d;
    long l;
    unsigned long ul;
    int8_t i8;
    int16_t i16;
    int32_t i32;
    int64_t i64;
    uint8_t ui8;
    uint16_t ui16;
    uint32_t ui32;
    uint64_t ui64;
    void *ptr;
  };
  Cache cache;
  
protected:
  Vertex *vertices[NVERT];
  Simplex *neighbors[NNEI]; 
  GlobalIdentity generation; // stores generation / index
  LocalIndex localIndex;
  Flag flags;

  void resetCache()
  {
    memset(cache.c,0,sizeof(cache));
  }
  
  Element* &getElementPtrRef(int id)
  {
    return vertices[id];
  }

  Element *getElementPtr(int id) const
  {
    return vertices[id];
  }

  void swapVertices(int a, int b)
  {
    std::swap(vertices[a],vertices[b]);	
    std::swap(neighbors[a],neighbors[b]);
  }

  void setSharedF(bool b=true)
  {
    if (b)
      flags |= SIMPLEX_FLAG_SHARED;
    else
      flags &= ~SIMPLEX_FLAG_SHARED;
  }

  void setBoundaryF(bool b=true)
  {
    if (b)
      flags |= SIMPLEX_FLAG_BOUNDARY;
    else
      flags &= ~SIMPLEX_FLAG_BOUNDARY;
  }
  
  void setSafeF(bool b=true)
  {
    if (b)
      flags &= ~SIMPLEX_FLAG_UNSAFE;
    else
      flags |= SIMPLEX_FLAG_UNSAFE;
  }

  void setSetF(bool b=true)
  {
    if (b)
      flags &= ~SIMPLEX_FLAG_NOTSET;
    else
      flags |= SIMPLEX_FLAG_NOTSET;
  }

  void setTaggedF(bool b=true)
  {
    if (b)
      flags |= SIMPLEX_FLAG_TAG;
    else
      flags &= ~SIMPLEX_FLAG_TAG;
  }

  
  template <class TT, class V>
  void setElements(TT* el[], const V &volumeFunctor)
  {
    Vertex *elements[NVERT];    
    for (int i=0;i<NVERT;i++) elements[i]=static_cast<Vertex*>(el[i]);
    setVertices(elements,volumeFunctor);    
  }

  template <class InputIterator, class V>
  void setVertices(InputIterator &vert,const V &volumeFunctor)
  {
    std::copy(vert,vert+NVERT,vertices);
    setSafeF(false);
    setSetF(true);
    orient(volumeFunctor); 
    setSafeF(true);
  }

  template <class InputIterator>
  void setVertices(InputIterator &vert)
  {
    std::copy(vert,vert+NVERT,vertices);
    setSafeF(false);
    setSetF(true);
  }
 
  void setVertex(int i, Vertex *v)
  {
    vertices[i]=v;
  }

  void setNeighbor(int i, Simplex *s)
  {
    neighbors[i]=s;
  }

  Simplex* checkNeighbors()
  {
    for (int i=0;i<NNEI;i++)
      {
	Simplex *s=neighbors[i];
	if (s!=NULL)
	  {
	    int count=0;
	    for (int j=0;j<NVERT;j++)
	      for (int k=0;k<NVERT;k++)
		{
		  if (vertices[j]==s->vertices[k])
		    count++;
		}
	    if (count != NVERT-1) return s;
	  }
      }
    return NULL;
  }

  // Ensures that the simplex is correctly oriented (i.e. volume is positive)  
  // orientation is determined by the NDIM first segments
  template <class V>
  bool orient(const V& volumeFunctor)
  {  
    if (isSafe()) return false;
  
    //Coord vec[NDIM][NDIM_W];
    //getBaseVectors(vec);
    
    double vol = volumeFunctor(static_cast<const Simplex *>(this));
    //volumeFunctor.get(static_cast<const Simplex *>(this));
    //printf("Vol = %lg\n",vol);swapVertices(NVERT-1,NVERT-2); 
    //printf("Vol2 = %lg\n",volumeFunctor(static_cast<const Simplex *>(this)));
    if (vol<0)
      swapVertices(NVERT-1,NVERT-2); 
    else 
      return false;
    
    return true;
  }
   
  void swap_internal(Simplex *other)
  {
    // PRESERVE THE SHADOW/GHOST FLAGS !!!!!!  
    Data::swap(other);
    std::swap(cache,other->cache);
    std::swap(vertices,other->vertices);
    std::swap(neighbors,other->neighbors);
    std::swap(localIndex,other->localIndex);  
    std::swap(*static_cast<Data*>(this),*static_cast<Data*>(other));
    Flag tmp(flags&(SIMPLEX_FLAG_SHADOW|SIMPLEX_FLAG_GHOST));
    flags=other->flags&(~(SIMPLEX_FLAG_SHADOW|SIMPLEX_FLAG_GHOST));   
    flags|=tmp;
  }
  
};

template <class T, template <class> class S>
const typename SimplexFromVerticesT<T,S>::LocalIndex
SimplexFromVerticesT<T,S>::LOCAL_INDEX_INVALID = T::LOCAL_INDEX_INVALID;

template <class T, template <class> class S>
const typename SimplexFromVerticesT<T,S>::LocalIndex
SimplexFromVerticesT<T,S>::LOCAL_INDEX_MAX = T::LOCAL_INDEX_MAX;

template <class T, template <class> class S>
const typename SimplexFromVerticesT<T,S>::GlobalIndex
SimplexFromVerticesT<T,S>::GLOBAL_INDEX_INVALID = T::GLOBAL_INDEX_INVALID;

template <class T, template <class> class S>
const typename SimplexFromVerticesT<T,S>::GlobalIndex
SimplexFromVerticesT<T,S>::GLOBAL_INDEX_MAX = T::GLOBAL_INDEX_MAX;

/** \}*/
#include "../internal/namespace.footer"
#endif
