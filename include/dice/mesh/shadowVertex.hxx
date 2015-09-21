#ifndef __SHADOW_VERTEX_HXX__
#define __SHADOW_VERTEX_HXX__

#include "../mesh/vertex.hxx"


#include "../internal/namespace.header"

template <class T> class LocalMeshT;
template <class T> class MeshT;
template <class T> class VertexT;

template <class T>
class ShadowVertexT : public VertexT<T>
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class VertexT<T>;
 
  typedef VertexT<T> Base;
  typedef ShadowVertexT<T> MyType;
  
  typedef typename Base::GlobalIndex GlobalIndex;  
  typedef typename Base::GlobalIdentity GlobalIdentity;
  typedef typename Base::GlobalIdentity::Rank GlobalIdentityRank;
  typedef typename Base::GlobalIdentity::Id GlobalIdentityId;
  
  ShadowVertexT():Base()
  {
    setShadowF(true); 
    Base::setSharedF(true); 
    Base::globalIdentity=GlobalIdentity::empty.get();
  }

  ~ShadowVertexT() {}

protected:

  void setShadowF(bool b=true)
  {
    if (b)
      Base::flags |= VERTEX_FLAG_SHADOW;
    else
      Base::flags &= ~VERTEX_FLAG_SHADOW;
  }
  
};

#include "../internal/namespace.footer"
#endif
