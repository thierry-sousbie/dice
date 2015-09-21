#ifndef __VERTEX_HXX__
#define __VERTEX_HXX__

#include <algorithm>

#include "../mesh/simplex.hxx"
#include "../mesh/vertexFlagsDefines.hxx"
#include "../mesh/ghostVertex.hxx"
#include "../mesh/shadowVertex.hxx"

#include "../tools/helpers/helpers.hxx"

/**
 * @file 
 * @brief Defines a vertex class to be used with MeshT
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

template <class T> class LocalMeshT;
template <class T> class MeshT;
template <class T> class SimplexT;

/**
 * \class VertexT
 * \brief Defines a vertex class to be used with MeshT
 * contrary to simplices, all Vertices store a local index AND a global identity 
 * because they may belong to the frontier of two nodes (i.e. local shared vertices).
 * As a consequence, the index part of the global identity of a regular vertex
 * cannot always be equal to its local index ...
 */
 

template <class T>
class VertexT : public T::VertexData
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class SimplexT<T>;
 
  static const long NDIM = T::NDIM;
  static const long NDIM_W = T::NDIM_W;
  
  typedef typename T::VertexFlag Flag;
  typedef typename T::VertexData Data;
  typedef typename T::Coord Coord;

  typedef typename T::LocalIndex LocalIndex;
  typedef typename T::GlobalIndex GlobalIndex;
  typedef typename T::GlobalIdentity GlobalIdentity;

  typedef VertexT<T> MyType;
  typedef T Traits; 
  
  typedef VertexT<T>       Vertex;
  typedef ShadowVertexT<T> Shadow;
  typedef GhostVertexT<T>  Ghost;
  
  VertexT():flags(VERTEX_FLAG_NOTSET)
  {
    //std::fill(coords,coords+NDIM_W,3);
    //std::fill(coords,coords+NDIM_W,0);
    //localIndex=0;
  }

  VertexT(const Coord *c, LocalIndex id):flags(VERTEX_FLAG_NOTSET)
  {
    set(c,id);
  }

  //static int getType() {return 0;}

  Coord *getCoordsPtr()
  {
    return coords;
  }  

  const Coord *getCoordsConstPtr() const
  {
    return coords;
  }  
  
  template <class OutputIterator, int D=NDIM_W>
  void getCoords(OutputIterator out)  const
  {
    for (int i=0;i<D;i++)
      {
	(*out)=coords[i];
	++out;
      }    
  }

  template <class InputIterator>
  void setCoords(InputIterator out)
  {
    for (int i=0;i<NDIM_W;i++)
      {
	coords[i]=(*out);
	++out;
      }    
  }

  template <class M, class SEG>
  static long getSimplexRefineBufferSize()
  {
    typedef typename M::Simplex Simplex;
    return T::template getSimplexRefineBufferSize<M,SEG,Vertex,Simplex>();
  }

  template <class M, class SEG>
  void setRefinedCoords(M *mesh, const SEG &seg, 
			typename M::Simplex * const *simplices, 
			int nSimplices,void* buffer)
  {
    return T::refineCoords(mesh,seg,this,simplices,nSimplices,buffer);
  }  
  
  Coord getCoord(int i)  const
  {
    return coords[i];
  }
  
  LocalIndex getLocalIndex() const
  {
    return localIndex;
  }
  
  GlobalIdentity getGlobalIdentity() const
  {
    return globalIdentity;
  }

  GlobalIdentity getGeneration() const
  {
    return generation;
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

  // void setData(const Data &src)
  // {
  //   return (*static_cast<const Data*>(this)) = src;
  // }

  bool isShadow() const
  {
    return (flags&VERTEX_FLAG_SHADOW);
  }

  bool isGhost() const
  {
    return (flags&VERTEX_FLAG_GHOST);
  }

  bool isShadowOrGhost() const
  {
    return flags & (VERTEX_FLAG_SHADOW|VERTEX_FLAG_GHOST);
  }

  bool isLocal() const
  {
    return !isShadowOrGhost();
  }

  bool isSet() const
  {
    return !(flags & VERTEX_FLAG_NOTSET);
  }

  bool isShared() const
  {
    return (flags & VERTEX_FLAG_SHARED);
  }

  bool isBoundary() const
  {
    return (flags & VERTEX_FLAG_BOUNDARY);
  }

  bool isTagged() const
  {
    return (flags&VERTEX_FLAG_TAG);
  }

  bool isTagged2() const
  {
    return (flags&VERTEX_FLAG_TAG2);
  }

  template <class L>
  void print() const
  {
    GlobalIdentity gid(getGlobalIdentity());
    long rank=gid.rank();
    long id=gid.id();
    if (NDIM_W==2)
      glb::console->print<L>("Vertex %ld(%ld) gid=(%ld,%ld) P=[%lg,%lg]:  Flags=%d\n",
			(long)this,(long)localIndex,rank,id,coords[0],coords[1],(int)flags);
//#if (DEFINED_MESH_DIM_WORLD>2)
    else if (NDIM_W==3)
      glb::console->print<L>("Vertex %ld(%ld) gid=(%ld,%ld) P=[%lg,%lg,%lg]: Flags=%d\n",
			(long)this,(long)localIndex,rank,id,coords[0],coords[1],coords[2],
			(int)flags);
//#if (DEFINED_MESH_DIM_WORLD>3)
    else if (NDIM_W==4)
      glb::console->print<L>("Vertex %ld(%ld) gid=(%ld,%ld) P=[%lg,%lg,%lg,%lg]: Flags=%d\n",
			(long)this,(long)localIndex,rank,id,coords[0],coords[1],coords[2],
			coords[3],(int)flags);
//#if (DEFINED_MESH_DIM_WORLD>4)
    else if (NDIM_W==5)
      glb::console->print<L>("Vertex %ld(%ld) gid=(%ld,%ld) P=[%lg,%lg,%lg,%lg,%lg]: Flags=%d\n",
			(long)this,(long)localIndex,rank,id,coords[0],coords[1],coords[2],
			coords[3],coords[4],(int)flags);
//#if (DEFINED_MESH_DIM_WORLD>5)
    else if (NDIM_W==6)
      glb::console->print<L>("Vertex %ld(%ld) gid=(%ld,%ld) P=[%.20lg,%.20lg,%.20lg,%lg,%lg,%lg]: Flags=%d\n",
			(long)this,(long)localIndex,rank,id,coords[0],coords[1],coords[2],
			coords[3],coords[4],coords[5],(int)flags);
//#endif //(DEFINED_MESH_DIM_WORLD>2)
//#endif //(DEFINED_MESH_DIM_WORLD>3)
//#endif //(DEFINED_MESH_DIM_WORLD>4)
//#endif //(DEFINED_MESH_DIM_WORLD>5)
  }

  Flag getFlags() const {return flags;}
  
  bool operator<  (const MyType &other) const {return globalIdentity< other.globalIdentity;}
  bool operator<= (const MyType &other) const {return globalIdentity<=other.globalIdentity;}
  bool operator>  (const MyType &other) const {return globalIdentity> other.globalIdentity;}
  bool operator>= (const MyType &other) const {return globalIdentity>=other.globalIdentity;}
  bool operator== (const MyType &other) const {return globalIdentity==other.globalIdentity;}
  bool operator!= (const MyType &other) const {return globalIdentity!=other.globalIdentity;}
  
  struct cmpPosLess
  {
    bool operator()(const MyType* a, const MyType *b) const
    {
      for (int i=0;i<NDIM_W;i++)
	{
	  if (a->coords[i]<b->coords[i]) return true;
	  else if (a->coords[i]>b->coords[i]) return false;
	}
      return false;
    }
  };

  template <class W>
  static void selfSerialize(const MyType *me, W *writer)
  {
    Data::selfSerialize(me,writer);
    writer->write(me->coords,NDIM_W);
    writer->write(&me->globalIdentity);
    writer->write(&me->generation);  
    writer->write(&me->localIndex);  
    writer->write(&me->flags);  
  }

  template <class R>
  static void selfUnSerialize(MyType *me, R *reader)
  {
    Data::selfUnSerialize(me,reader);
    reader->read(me->coords,NDIM_W);
    reader->read(&me->globalIdentity);
    reader->read(&me->generation);  
    reader->read(&me->localIndex);  
    reader->read(&me->flags);  
  }

  // union Cache
  // {
  //   struct {float f; int i;} pfi;
  //   int i[2];
  //   float f[2];
  //   double d;
  //   long l;
  //   unsigned long ul;
  //   void *ptr;
  // };
  // Cache cache;

protected:  
  Coord coords[NDIM_W];
  GlobalIdentity globalIdentity;
  GlobalIdentity generation; // stores generation / index
  LocalIndex localIndex;
  Flag flags;

  // this will conserve the shadow/ghost status of the vertex
  void copy(const Vertex *v)
  {
    Flag tmp(flags&(VERTEX_FLAG_SHADOW|VERTEX_FLAG_GHOST));
    if (tmp) tmp|=VERTEX_FLAG_SHARED; // shadows and ghosts are always shared !
    std::copy(v->coords,v->coords+NDIM_W,coords);
    setData(v->getData());//*static_cast<Data*>(this) = *static_cast<Data*>(v) ;  
    globalIdentity = v->globalIdentity();
    localIndex=v->localIndex;
    generation=v->generation;
    flags=v->flags&(~(VERTEX_FLAG_SHADOW|VERTEX_FLAG_GHOST));   
    flags|=tmp;   
  }
  
  void setLocalIndex(LocalIndex i)
  {
    localIndex = i;
  }

  void setGeneration(long gen,unsigned long index)
  {    
    generation = GlobalIdentity(gen,index);
  }

  void setGeneration(GlobalIdentity gen)
  {    
    generation = gen;
  }

  void setGlobalIdentity(long node, long index)
  {    
    globalIdentity = GlobalIdentity(node,index);
  }

  void setGlobalIdentity(GlobalIdentity id_)
  {
    globalIdentity = id_;
  }

  void set(const Coord *c, LocalIndex id)
  {
    std::copy(c,c+NDIM_W,coords);
    localIndex=id; 
    setSetF(true);
  }

  void setSetF(bool b=true)
  {
    if (b)
      flags &= ~VERTEX_FLAG_NOTSET;
    else
      flags |= VERTEX_FLAG_NOTSET;
  }

  void setSharedF(bool b=true)
  {
    if (b)
      flags |= VERTEX_FLAG_SHARED;
    else
      flags &= ~VERTEX_FLAG_SHARED;
  }
 
  void setBoundaryF(bool b=true)
  {
    if (b)
      flags |= VERTEX_FLAG_BOUNDARY;
    else
      flags &= ~VERTEX_FLAG_BOUNDARY;
  }

  void setTaggedF(bool b=true)
  {
    if (b)
      flags |= VERTEX_FLAG_TAG;
    else
      flags &= ~VERTEX_FLAG_TAG;
  }

  void setTagged2F(bool b=true)
  {
    if (b)
      flags |= VERTEX_FLAG_TAG2;
    else
      flags &= ~VERTEX_FLAG_TAG2;
  }

  void cleanTags()
  {
    flags &=(~(VERTEX_FLAG_TAG|VERTEX_FLAG_TAG2));
  }

};

/** \}*/
#include "../internal/namespace.footer"
#endif
