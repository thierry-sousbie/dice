#ifndef __GHOST_VERTEX_HXX__
#define __GHOST_VERTEX_HXX__

#include "../mesh/vertex.hxx"


#include "../internal/namespace.header"

template <class T> class LocalMeshT;
template <class T> class MeshT;
template <class T> class VertexT;

template <class T>
class GhostVertexT : public VertexT<T>
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class VertexT<T>;
 
  typedef VertexT<T> Base;
  typedef GhostVertexT<T> MyType;
  
  typedef typename Base::GlobalIndex GlobalIndex;  
  typedef typename Base::GlobalIdentity GlobalIdentity;
  typedef typename Base::GlobalIdentity::Rank GlobalIdentityRank;
  typedef typename Base::GlobalIdentity::Id GlobalIdentityId;
  
  GhostVertexT():Base()
  {
    setGhostF(true); 
    Base::setSharedF(true); 
    Base::globalIdentity=GlobalIdentity::empty.get();
  }

  ~GhostVertexT() {}

protected:

  void setGhostF(bool b=true)
  {
    if (b)
      Base::flags |= VERTEX_FLAG_GHOST;
    else
      Base::flags &= ~VERTEX_FLAG_GHOST;
  }
  
};

#include "../internal/namespace.footer"
#endif
