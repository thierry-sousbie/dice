#ifndef __SHADOW_SIMPLEX_HXX__
#define __SHADOW_SIMPLEX_HXX__

#include "../mesh/simplex.hxx"

/**
 * @file 
 * @brief Defines shadow simplices class.
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
 * \class ShadowSimplexT
 * \brief A shadow simplex is a simplex that is the image of a simplex local to a remote
 * MPI node and that does not contain any local vertex. Contrary to ghost simplices 
 * (see GhostSimplexT), only Shadow simplices neighbors of the global mesh that 
 * are also Ghost simplices are guaranteed to be locally accessible.
 */ 

template <class T>
class ShadowSimplexT : public SimplexT<T>
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class SimplexT<T>;
 
  typedef SimplexT<T> Base;
  typedef typename SimplexT<T>::LeafBase LeafBase;
  typedef ShadowSimplexT<T> MyType;

  typedef typename Base::Segment Segment;
  typedef typename Base::Vertex Vertex; 
  typedef typename Base::GlobalIndex GlobalIndex;  
  typedef typename Base::Element Element;
  typedef typename Base::GlobalIdentity GlobalIdentity;
  typedef typename Base::GlobalIdentity::Rank GlobalIdentityRank;
  typedef typename Base::GlobalIdentity::Id GlobalIdentityId;

  ShadowSimplexT():Base()
  {
    setShadowF(true); 
    Base::setSharedF(true);
    globalIdentity=GlobalIdentity::empty.get();
  }

  ~ShadowSimplexT() {}

  template <class W>
  static void selfSerialize(const MyType *me, W *writer)
  {
    writer->write(&me->globalIdentity);
    Base::selfSerialize(static_cast<const Base*>(me), writer);
  }

  template <class R>
  static void selfUnSerialize(MyType *me, R *reader)
  {
    reader->read(&me->globalIdentity);
    Base::selfUnSerialize(static_cast<Base*>(me), reader);
  }

  GlobalIdentity getGlobalIdentity() const
  {   
    return globalIdentity;
  }

  template <class L>
  void print()
  {
    char tmp[60];
    //sprintf(tmp,"Shadow(%ld,%ld)",(long)globalIdentity.rank(),(long)globalIdentity.id());
    sprintf(tmp,"Shadow@");
    Base::template print<L>(std::string(tmp));
  }

  friend struct cmpPtrMore;
  friend struct cmpPtrLess;
  struct cmpPtrLess
  {
    bool operator()(const MyType* a,const MyType* b)
    {
      return a->globalIdentity<b->globalIdentity;
    }
  };
  struct cmpPtrMore
  {
    bool operator()(const MyType* a,const MyType* b)
    {
      return a->globalIdentity>b->globalIdentity;
    }
  };
  /*
  static bool cmpGidLess(const MyType& a, const MyType& b)
  {
    return a.globalIdentity<b.globalIdentity;
  }

  static bool cmpGidMore(const MyType& a, const MyType& b)
  {
    return a.globalIdentity>b.globalIdentity;
    }
  */

  /*
  bool operator<  (const MyType &other) const {return globalIdentity<other.globalIdentity;}
  bool operator<= (const MyType &other) const {return globalIdentity<=other.globalIdentity;}
  bool operator>  (const MyType &other) const {return globalIdentity>other.globalIdentity;}
  bool operator>= (const MyType &other) const {return globalIdentity>=other.globalIdentity;}
  bool operator== (const MyType &other) const {return globalIdentity==other.globalIdentity;}
  bool operator!= (const MyType &other) const {return globalIdentity!=other.globalIdentity;}
  */
protected:
  GlobalIdentity globalIdentity;

  // Shadows have no tree structure, but they still have a parent pointer, so we can
  // use it internally as a temporary variable, to store a pointer to a partner simplex
  // void setTemporaryPartner(MyType *simplex)
  // {
  //   LeafBase::parent = static_cast<LeafBase*>(simplex);
  // }

  template <class TT, class V>
  void setElements(TT* el[], const V &volumeFunctor)
  {
    Base::setElements(el,volumeFunctor);
  }

  void setGlobalIdentity(GlobalIdentityRank n, GlobalIdentityId id)
  {
    globalIdentity.set(n,id);
  }

  void setGlobalIdentity(GlobalIdentity id)
  {
    globalIdentity=id;
  }
  
  void setShadowF(bool b=true)
  {
    if (b)
      Base::flags |= SIMPLEX_FLAG_SHADOW;
    else
      Base::flags &= ~SIMPLEX_FLAG_SHADOW;
  }
  
  void addNeighbor(Vertex *vertex, Base *simplex)
  {
    int i;
    for (i=0;i<Base::NVERT;i++)
      if (Base::getElementPtr(i)==static_cast<Element*>(vertex)) break;
    Base::neighbors[i]=simplex;
  }

  void setNeighbor(int index, Base *simplex)
  {
    Base::neighbors[index]=simplex;
  }

  template <class VU, class GVU, class SVU, class SU, class GSU, class SSU>
  void updateAfterUnserialized(const VU  &vertexUpdate,
			       const GVU &ghostVertexUpdate,
			       const SVU &shadowVertexUpdate,
			       const SU  &simplexUpdate,
			       const GSU &ghostSimplexUpdate,
			       const SSU &shadowSimplexUpdate,
			       bool swap=false)
  {
    Base::updateAfterUnserialized(vertexUpdate,
				  ghostVertexUpdate, shadowVertexUpdate,
				  simplexUpdate,
				  ghostSimplexUpdate,shadowSimplexUpdate,
				  swap);
  }
  
};

/** \}*/
#include "../internal/namespace.footer"
#endif
