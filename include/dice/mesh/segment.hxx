#ifndef __SEGMENT_HXX__
#define __SEGMENT_HXX__

#include <algorithm>

#include "../mesh/vertex.hxx"
#include "../mesh/simplex.hxx"

/**
 * @file
 * @brief Defines a segment class to be used with MeshT
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

template <class T>
class VertexT;
template <class T>
class SimplexT;
template <class T>
class LocalMeshT;
template <class T>
class MeshT;

#define SEGMENT_FLAG_NOTSET (1 << 0)
#define SEGMENT_FLAG_ALERT (1 << 7)
// #define SEGMENT_FLAG_SHARED   (1<<1)
// #define SEGMENT_FLAG_BOUNDARY (1<<2)

/**
 * \class SegmentT
 * \brief Defines a segment class to be used with MeshT.
 */

template <class T>
class SegmentT
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class SimplexT<T>;
  friend class VertexT<T>;

  static const long NDIM = T::NDIM;
  static const long NDIM_W = T::NDIM_W;
  static const long NVERT = 2;

  typedef typename T::SegmentFlag Flag;

  typedef SegmentT<T> MyType;
  typedef T Traits;

  typedef SimplexT<T> Simplex;
  typedef SegmentT<T> Segment;
  typedef VertexT<T> Vertex;

  typedef typename Vertex::GlobalIdentity GlobalIdentity;
  typedef typename Vertex::LocalIndex LocalIndex;

  typedef typename Simplex::SegmentHandle Handle;
  typedef typename Simplex::ConstSegmentHandle ConstHandle;
  typedef typename Simplex::FacetHandle FacetHandle;
  typedef typename Simplex::ConstFacetHandle ConstFacetHandle;

  typedef typename Simplex::segment_circulator circulator;
  typedef typename Simplex::const_segment_circulator const_circulator;

  SegmentT() : flags(SEGMENT_FLAG_NOTSET)
  {
    vertices[0] = vertices[1] = NULL;
    simplex = NULL;
  }

  SegmentT(Vertex *v[2], Simplex *s = NULL) : flags(SEGMENT_FLAG_NOTSET)
  {
    set(v, s);
  }

  SegmentT(Vertex *v0, Vertex *v1, Simplex *s = NULL) : flags(SEGMENT_FLAG_NOTSET)
  {
    set(v0, v1, s);
  }

  long getGeneration() const
  {
    return 1 + std::max(vertices[0]->getGeneration().rank(),
                        vertices[1]->getGeneration().rank());
  }

  // static int getType() {return 1;}

  Simplex *getSimplex() const
  {
    return simplex;
  }

  Vertex *getOtherVertex(Vertex *v) const
  {
    if (v == vertices[0])
      return vertices[1];
    if (v == vertices[1])
      return vertices[0];
    return NULL;
  }

  Vertex *getLowVertex() const
  {
    if ((*vertices[0]) < (*vertices[1]))
      return vertices[0];
    else
      return vertices[1];
  }

  Vertex *getHighVertex() const
  {
    if ((*vertices[1]) < (*vertices[0]))
      return vertices[0];
    else
      return vertices[1];
  }

  std::pair<Vertex *, Vertex *> getOrderedVertices() const
  {
    if ((*vertices[0]) < (*vertices[1]))
      return std::make_pair(vertices[0], vertices[1]);
    else
      return std::make_pair(vertices[1], vertices[0]);
  }

  Vertex *getVertex(int which, bool reverse) const
  {
    return (reverse) ? vertices[1 - which] : vertices[which];
  }

  Vertex *getVertex(int which) const
  {
    return vertices[which];
  }

  Vertex *getVertexByGlobalId(GlobalIdentity vid) const
  {
    if (vertices[0]->getGlobalIdentity() == vid)
      return vertices[0];
    if (vertices[1]->getGlobalIdentity() == vid)
      return vertices[1];
    return NULL;
  }

  Vertex *getVertexByLocalIndex(LocalIndex vid) const
  {
    if (vertices[0]->getLocalIndex() == vid)
      return vertices[0];
    if (vertices[1]->getLocalIndex() == vid)
      return vertices[1];
    return NULL;
  }

  template <class OutputIterator>
  void getVertices(OutputIterator out, bool reverse) const
  {
    if (reverse)
    {
      *out = vertices[1];
      ++out;
      *out = vertices[0];
    }
    else
    {
      *out = vertices[0];
      ++out;
      *out = vertices[1];
    }
  }

  template <class OutputIterator>
  void getVertices(OutputIterator out) const
  {
    *out = vertices[0];
    ++out;
    *out = vertices[1];
  }

  template <class OutputIterator>
  void getVerticesLocalIndex(OutputIterator out, bool reverse) const
  {
    if (reverse)
    {
      *out = vertices[1]->getLocalIndex();
      ++out;
      *out = vertices[0]->getLocalIndex();
    }
    else
    {
      *out = vertices[0]->getLocalIndex();
      ++out;
      *out = vertices[1]->getLocalIndex();
    }
  }

  template <class OutputIterator>
  void getVerticesLocalIndex(OutputIterator out) const
  {
    *out = vertices[0]->getLocalIndex();
    ++out;
    *out = vertices[1]->getLocalIndex();
  }

  circulator getCirculator() const
  {
    // Handle handle(this);
    return circulator(Handle(this), simplex);
  }

  circulator getCirculator(Vertex *v) const
  {
    // Handle handle(this);
    return circulator(Handle(this), v, simplex);
  }

  const_circulator getConstCirculator() const
  {
    // Handle handle(this);
    return const_circulator(Handle(this), simplex);
  }

  const_circulator getConstCirculator(Vertex *v) const
  {
    // Handle handle(this);
    return const_circulator(Handle(this), v, simplex);
  }

  template <class OutputIterator>
  int getAdjacentSimplices(OutputIterator out)
  {
    circulator ci_end = circulator(Handle(this), simplex);
    circulator ci = ci_end;
    int j = 0;
    do
    {
      *out = *ci;
      ++j;
      ++out;
    } while ((++ci) != ci_end);

    return j;
  }

  bool isSet() const
  {
    return !(flags & SEGMENT_FLAG_NOTSET);
  }

  bool isShared() const
  {
    return (vertices[0]->isShared() || vertices[1]->isShared());
  }

  bool isOnBoundary() const
  {
    return (vertices[0]->isBoundary() && vertices[1]->isBoundary());
  }

  bool isGhost() const
  {
    return (vertices[0]->isGhost() || vertices[1]->isGhost());
  }

  bool isShadow() const
  {
    return (vertices[0]->isShadow() || vertices[1]->isShadow());
  }

  bool isShadowOrGhost() const
  {
    return (vertices[0]->isShadowOrGhost() || vertices[1]->isShadowOrGhost());
  }

  bool haveAlertFlag() const
  {
    return (flags & SEGMENT_FLAG_ALERT);
  }

  /*
  bool isShared() const
  {
    return flags & SEGMENT_FLAG_SHARED;
  }

  bool isBoundary() const
  {
    return flags & SEGMENT_FLAG_BOUNDARY;
  }
  */

  /*
  // FIXME boundary conditons
  double length() const
  {
    double len=0;
    for (int i=0;i<NDIM_W;i++)
      {
  double d=(vertices[1]->coords[i]-vertices[0]->coords[i]);
  if (d>=1) d-=2;
  if (d<-1) d+=2;
  len+=d*d;
      }
    return sqrt(len);
  }
  */
  Handle getHandle() const
  {
    return Handle(this);
  }

  ConstHandle getConstHandle() const
  {
    return ConstHandle(this);
  }

  int countSimplices() const
  {
    if (NDIM < 2)
      return 0;
    // Handle handle(this);
    int N = 0;
    // simplex_circulator ci_end=simplex->getSegmentCirculator(handle);
    const_circulator ci_end = getCirculator();
    const_circulator ci = ci_end;

    do
    {
      N++;
    } while ((++ci) != ci_end);

    return N;
  }

  // returns true if there was any non local simplex
  template <class IT>
  bool countLocalSimplices(IT &val) const
  {
    if (NDIM < 2)
    {
      val = 0;
      return false;
    }
    // Handle handle(this);
    // simplex_circulator ci_end=simplex->getSegmentCirculator(handle);
    const_circulator ci_end = getConstCirculator();
    const_circulator ci = ci_end;
    bool found = false;

    val = 0;
    do
    {
      if (ci->isShadowOrGhost())
        found = true;
      else
        ++val;
    } while ((++ci) != ci_end);

    return found;
  }

  /*
  template <class OutputIterator>
  void getCenter(OutputIterator c) const
  {
    for (int i=0;i<NDIM_W;i++)
      {
  (*c)=0.5*(vertices[1]->getCoord(i)+vertices[0]->getCoord(i));
  ++c;
      }
  }
  */
  template <class L>
  void print() const
  {
    glb::console->print<L>("Segment [%ld(%ld),%ld(%ld)]=[(%ld,%ld);(%ld,%ld)]: s=%ld(%ld) / Flags=%d\n",
                           (long)vertices[0], (long)vertices[0]->getLocalIndex(),
                           (long)vertices[1], (long)vertices[1]->getLocalIndex(),
                           (long)vertices[0]->getGlobalIdentity().rank(), (long)vertices[0]->getGlobalIdentity().id(),
                           (long)vertices[1]->getGlobalIdentity().rank(), (long)vertices[1]->getGlobalIdentity().id(),
                           (long)simplex, (long)simplex->getLocalIndex(), (int)flags);
  }

  bool operator<(const MyType &other) const
  {
    std::pair<Vertex *, Vertex *> A = getOrderedVertices();
    std::pair<Vertex *, Vertex *> B = other.getOrderedVertices();
    if ((*A.first) < (*B.first))
      return true;
    else if ((*A.first) == (*B.first))
    {
      if ((*A.second) < (*B.second))
        return true;
    }

    return false;
  }

protected:
  Vertex *vertices[2];
  Simplex *simplex;
  Flag flags;

  void split(Vertex *v, Segment *newSeg)
  {
    (*newSeg) = (*this);
    vertices[1] = v;
    newSeg->vertices[0] = v;
  }

  void setVertex(int which, Vertex *v, bool reverse = false)
  {
    vertices[(reverse) ? (1 - which) : which] = v;
  }

  void set(Vertex *v[2], Simplex *s = NULL)
  {
    vertices[0] = v[0];
    vertices[1] = v[1];
    setSimplex(s);
  }

  void set(Vertex *v0, Vertex *v1, Simplex *s = NULL)
  {
    vertices[0] = v0;
    vertices[1] = v1;
    setSimplex(s);
  }

  void setSimplex(Simplex *s)
  {
    simplex = s;
    if (s != NULL)
    {
      if ((vertices[0] != NULL) && (vertices[1] != NULL))
        setSetF(true);
    }
  }

  void setSetF(bool b = true)
  {
    if (b)
      flags &= ~SEGMENT_FLAG_NOTSET;
    else
      flags |= SEGMENT_FLAG_NOTSET;
  }

  void setAlertF(bool b = true)
  {
    if (b)
      flags |= SEGMENT_FLAG_ALERT;
    else
      flags &= ~SEGMENT_FLAG_ALERT;
  }
  /*
  void setBoundaryF(bool b=true)
  {
    if (b)
      flags |= SEGMENT_FLAG_BOUNDARY;
    else
      flags &= ~SEGMENT_FLAG_BOUNDARY;

    if (vertices[0]!=NULL) vertices[0]->setBoundaryF(b);
    if (vertices[1]!=NULL) vertices[1]->setBoundaryF(b);
  }

  void setSharedF(bool b=true)
  {
    if (b)
      flags |= SEGMENT_FLAG_SHARED;
    else
      flags &= ~SEGMENT_FLAG_SHARED;

    if (vertices[0]!=NULL) vertices[0]->setSharedF(b);
    if (vertices[1]!=NULL) vertices[1]->setSharedF(b);
  }
  */
};

/** \}*/
#include "../internal/namespace.footer"
#endif
