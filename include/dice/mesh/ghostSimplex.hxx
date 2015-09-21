#ifndef __GHOST_SIMPLEX_HXX__
#define __GHOST_SIMPLEX_HXX__

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
 * \class GhostSimplexT
 * \brief A ghost simplex is a simplex that is the image of a simplex local to a remote
 * MPI node. While any remote simplex may be stored as a ghost simplex, any simplex 
 * containing a vertex local to current MPI node must be a ghost simplex. Contrary to
 * shadow simplices (see ShadowSimplexT), all ghost simplices neighbors existing in the 
 * global mesh are guaranteed to be also be stored locally as local, ghost or shadow 
 * simplices.
 */ 


template <class T>
class GhostSimplexT : public SimplexT<T>
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class SimplexT<T>;
 
  typedef SimplexT<T> Base;
  typedef typename SimplexT<T>::LeafBase LeafBase;
  typedef GhostSimplexT<T> MyType;

  typedef typename Base::Segment Segment;
  typedef typename Base::Vertex Vertex; 
  typedef typename Base::GlobalIndex GlobalIndex;  
  typedef typename Base::Element Element;
  typedef typename Base::GlobalIdentity GlobalIdentity;
  typedef typename Base::GlobalIdentity::Rank GlobalIdentityRank;
  typedef typename Base::GlobalIdentity::Id GlobalIdentityId;

  GhostSimplexT():Base()
  {
    setGhostF(true); 
    Base::setSharedF(true);
    globalIdentity=GlobalIdentity::empty.get();
  }

  ~GhostSimplexT() {}

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
    //sprintf(tmp,"Ghost(%ld,%ld)",(long)globalIdentity.rank(),(long)globalIdentity.id());
    sprintf(tmp,"Ghost@");
    Base::template print<L>(std::string(tmp));
  } 

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

  // Ghosts have no tree structure, but they still have a parent pointer, so we can
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
  
  void setGhostF(bool b=true)
  {
    if (b)
      Base::flags |= SIMPLEX_FLAG_GHOST;
    else
      Base::flags &= ~SIMPLEX_FLAG_GHOST;
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
