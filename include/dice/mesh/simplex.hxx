#ifndef __SIMPLEX_HXX__
#define __SIMPLEX_HXX__

#include "../mesh/internal/segmentCirculator.hxx"

#include "../mesh/simplexType.hxx"
#include "../mesh/sharedTree.hxx"
#include "../mesh/sharedTreeNodes.hxx"
#include "../mesh/simplexFlagsDefines.hxx"
#include "../mesh/simplexFromVertices.hxx"
#include "../mesh/ghostSimplex.hxx"
#include "../mesh/shadowSimplex.hxx"

#include "../geometry/predicates/pointInSimplex.hxx"

/**
 * @file 
 * @brief Defines Simplex data type.
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
 * \class SimplexBaseT
 * \brief base for SimplexT class defining the actual implementation of the simplex. This
 * class can be specialized to a particular implement using the second template parameter.
 * \note At the moment, only simplexType::VerticesOnly works, that is SimplexBaseT
 * must be specialized to SimplexFromVerticesT.
 * \tparam T mesh traits. See MeshTraitsT.
 * \tparam TT specialization. At the moment, only simplexType::VerticesOnly works, which
 * is equivalent to unsing SimplexFromVerticesT.
 */ 
template <class T, int TT> class SimplexBaseT;

template <class T> 
class SimplexBaseT< T , simplexType::VerticesOnly > : 
  public SimplexFromVerticesT<T,SimplexT> 
{
};

/**
 * \class SimplexT 
 * \brief Simplex data type to be used as cells in MeshT. This is just an interface, see
 * also SimplexBaseT.
 * \tparam mesh traits. See MeshTraitsT.
 */

template <class T>
class SimplexT : 
  public sharedTreeNodes::LeafT<SimplexT<T>,T>,
//public SimplexFromVerticesT< T, SimplexT >
public SimplexBaseT< T , T::SIMPLEX_TYPE >
  
{
public:
  friend class LocalMeshT<T>;
  friend class MeshT<T>;
  friend class VertexT<T>; 
  
  typedef T Traits;
  typedef sharedTreeNodes::LeafT<SimplexT<T>,T> LeafBase;
  typedef SimplexBaseT<T,T::SIMPLEX_TYPE>  Base;
  //typedef SimplexFromVerticesT<T,SimplexT> Base;

  static const int NDIM = Base::NDIM;
  static const int NDIM_W = Base::NDIM_W;
  static const int NNEI = Base::NNEI; 
  static const int NVERT = Base::NVERT;
  static const int NSEG = Base::NSEG;
  static const int NFACET = Base::NFACET;
  static const int NOPPEL = Base::NOPPEL; 

  static const bool EXPLICIT_SEGMENTS = Base::EXPLICIT_SEGMENTS;  
  
  typedef SimplexT<T>       MyType;
  typedef SimplexT<T>       Simplex;
  typedef ShadowSimplexT<T> Shadow;
  typedef GhostSimplexT<T>  Ghost;
 
  typedef typename Base::Vertex        Vertex;
  typedef typename Base::Element       Element;
  typedef typename Base::Segment       Segment;
  typedef typename Base::Facet         Facet;
  typedef typename Base::SegmentHandle SegmentHandle;
  typedef typename Base::FacetHandle   FacetHandle;

  typedef typename Vertex::Shadow      ShadowVertex;
  typedef typename Vertex::Ghost       GhostVertex;

  typedef internal::SegmentCirculatorT<SegmentHandle,Simplex> segment_circulator;
  typedef internal::SegmentCirculatorT<const SegmentHandle,const Simplex> const_segment_circulator;
  
  typedef typename T::LocalIndex LocalIndex;
  typedef typename T::GlobalIndex GlobalIndex;
  typedef typename T::GlobalIdentity GlobalIdentity;
  typedef typename T::Coord Coord;

  struct Neighborhood {
    typedef MyType Simplex;
    typedef typename Simplex::Vertex Vertex;
    typedef typename Simplex::SegmentHandle SegmentHandle;
    typedef typename Simplex::FacetHandle FacetHandle;

    FacetHandle fh[Simplex::NFACET];   /**< handles to NDIM+1 facets */
    Simplex *neiS[Simplex::NNEI];      /**< pointers to neighbor simplices (same order as facets)*/
    Vertex  *curV[Simplex::NVERT];     /**< vertices opposite to each facet in the current simplex (same order as facets)*/
    Vertex  *neiV[Simplex::NNEI];      /**< vertices opposite to each facet in the neighbor simplices (one for each neighbor, in the same order as facets)*/
    Simplex *simplex;     /**< the current simplex from which the neighborhood is computed */
    char nNeighbors;      /**< the number of neighbors (<=Simplex::NNEI, depending on the boundary conditions)*/

    Neighborhood()
    {nNeighbors=0;}
    explicit Neighborhood(Simplex *s)    
    {set(s);}  
    void set(Simplex *s)
    {s->getNeighborhood(*this);}
  };

  SimplexT():
    LeafBase(),
    Base() 
  {

  }

  static int getType() {return Base::NDIM;}

  template <class W>
  static void selfSerialize(const MyType *me, W *writer)
  {
    LeafBase::selfSerialize(static_cast<const LeafBase*>(me), writer);
    Base::selfSerialize(static_cast<const Base*>(me), writer);
  }

  template <class R>
  static void selfUnSerialize(MyType *me, R *reader)
  {
    LeafBase::selfUnSerialize(static_cast<LeafBase*>(me), reader);
    Base::selfUnSerialize(static_cast<Base*>(me), reader);
  }

  const Shadow *getAsShadow() const
  {
    if (!Base::isShadow()) 
      {
	//glb::console->print<LOG_WARNING>("calling getAsShadow() on a non shadow simplex!\n");
	return NULL;
      }
    else 
      return static_cast< const Shadow* >(this);
  }

  const Ghost *getAsGhost() const
  {
    if (!Base::isGhost()) 
      return NULL;
    else 
      return static_cast< const Ghost* >(this);
  }

  Shadow *getAsShadow()
  {
    if (!Base::isShadow()) 
      return NULL;
    else 
      return static_cast< Shadow* >(this);
  }

  Ghost *getAsGhost()
  {
    if (!Base::isGhost()) 
      return NULL;
    else 
      return static_cast< Ghost* >(this);
  }

  segment_circulator getSegmentCirculator(SegmentHandle &seg)
  {
    return segment_circulator(seg,this);
  }
  
  const_segment_circulator getSegmentCirculator(const SegmentHandle &seg) const
  {
    return const_segment_circulator(seg,this);
  }

  int getNeighborhood(Neighborhood &nb)
  {
    nb.nNeighbors = 0;
    nb.simplex = this;
    for (int i=0;i<NFACET;++i)
      {
	nb.curV[i] = this->getVertex(i);
	nb.fh[i] = this->getFacetHandle(i);
	if (nb.fh[i]->swapSimplex() != NULL)
	  {
	    nb.neiS[i]=nb.fh[i]->getSimplex();
	    nb.neiV[i]=nb.fh[i]->getOppositeVertex();
	    nb.nNeighbors++;
	  }
	else 
	  {
	    nb.neiV[i]=NULL;
	    nb.neiS[i]=NULL;
	  }
      }

    return nb.nNeighbors;
  }
  
  template <typename OutputIterator>
  int getOwnedLocalFacetHandles(OutputIterator handles, bool noGhost)
  {
    int count=0;
    Simplex *nei[NNEI];
    
    if (noGhost)
      {
	if (!Base::isShadowOrGhost())
	  {
	    for (int i=0;i<NFACET;++i)
	      {
		if ((Base::neighbors[i]==NULL)||
		    (Base::neighbors[i]->isShadowOrGhost())||
		    (Base::neighbors[i] > this)
		    //(Base::neighbors[i]->getLocalIndex() > getLocalIndex())
		    )
		  {
		    (*handles) = Base::getFacetHandle(i);
		    ++handles;
		    ++count;
		  }
	      }
	  }
      }
    else
      {
	if (!Base::isShadow())
	  {
	    for (int i=0;i<NFACET;++i)
	      {
		if ((Base::neighbors[i]==NULL)||
		    (Base::neighbors[i]->isShadow())||
		    (Base::neighbors[i] > this)
		    //(Base::neighbors[i]->getGlobalIndex() > getGlobalIndex())
		    )
		  {
		    (*handles) = Base::getFacetHandle(i);
		    ++handles;
		    ++count;
		  }
	      }
	  }
      }

    return count;
  }

  template <typename OutputIterator>
  int getOwnedLocalSegmentHandles(OutputIterator handles, bool noGhost)
  {
    char tag;
    return getOwnedLocalSegmentHandles(handles,noGhost,tag);
  }
  
  template <typename OutputIterator>
  int getOwnedLocalSegmentHandles(OutputIterator handles, bool noGhost, char &tag)
  {
    int count=0;
    tag=0;

    if (noGhost)
      {
	// Shadow and ghost simplices do not own any segment so that we only 
	// need to iterate over the local simplices
	if (!Base::isShadowOrGhost())
	  {
	    for (int i=0;i<NSEG;++i)
	      {
		SegmentHandle sh=Base::getSegmentHandle(i);
		segment_circulator sc = sh->getCirculator();	
		do {
		  ++sc;
		} while (((*sc) < this)||(sc->isShadowOrGhost()));
		
		// If no other local simplex containing segment 'i' has a smaller pointer 
		// => we own it
		if ((*sc)==this) 
		  {
		    (*handles) = sh;
		    ++handles;
		    ++count;
		    tag|=(1<<i);
		  }
	      }
	  }
      }
    else
      {
	// Shadow simplices do not own any segment 
	if (!Base::isShadow())
	  {
	    for (int i=0;i<NSEG;++i)
	      {
		SegmentHandle sh=Base::getSegmentHandle(i);
		segment_circulator sc = sh->getCirculator();	
		do {
		  ++sc;
		} while (((*sc) < this)||(sc->isShadow()));
		
		// If no other local/ghost simplex containing segment 'i' has a 
		// smaller pointer => we own it	
		if ((*sc)==this) 
		  {
		    (*handles) = sh;
		    ++handles;
		    ++count;
		    tag|=(1<<i);
		  }
	      }
	  }
      }

    return count;
  } 

  template <typename OutputIterator>
  int getOwnedLocalSegmentHandlesUseTag(OutputIterator handles, const char tag)
  {
    int count=0;

    if (tag>0)
      {
	for (int i=0;i<NSEG;++i)
	  {
	    if (tag&(1<<i))
	      {
		(*handles) = Base::getSegmentHandle(i);
		++handles;
		++count;
	      }
	  }
      }

    return count;
  }
  /*
  template <typename OutputIterator>
  int getVertexStar(int vertexId, OutputIterator out) const
  {
    int count=0;

    for (int i=0;i<NNEI;++i)
      {
	if (i!=vertexId) 
	  {
	    if (Base::neighbors[i]!=NULL)
	      {
		(*out) = Base::neighbors[i];
		++out;
		++count;
		Base::neighbors[i]->getVertexStar(vertexId,out);
		
	      }
	  }
      }

    return count;
  }
  */
  Simplex *getNeighbor(int i) const
  {
    // if ((i<0)||(i>NNEI))
    //   {
    // 	glb::console->print<LOG_ERROR>("Index out of range : %d\n",i);
    // 	exit(0);
    //   }
    return Base::neighbors[i];
  }

  // get all the vertices except 'ref'
  template <typename OutputIterator>
  int getNeighbors(OutputIterator out) const
  {
    int count=0;
    for (int i=0;i<NNEI;i++)
      {
	if (Base::neighbors[i]==NULL) continue;
	(*out)=Base::neighbors[i];
	++out;
	++count;
      }
    return count;
  }

  const Simplex** getNeighborsArr() const
  {
    return Base::neighbors;
  }

  // get all the vertices except 'ref'
  template <typename OutputIterator>
  int getNeighborsLocalIndex(OutputIterator out) const
  {
    int count=0;
    for (int i=0;i<NNEI;i++)
      {
	if (Base::neighbors[i]==NULL) continue;
	(*out)=Base::neighbors[i]->getLocalIndex();
	++out;
	++count;
      }
    return count;
  }
  /*
  // return +1 if same (or does not exist), -1 if not
  int compareNeighborOrientation(int index)
  {
    Simplex *nei = getNeighbor(index);
    if (nei==NULL) return 1;
    int oppIndex = nei->getNeighborIndex(this);
    int order[NVERT];

    for (int i=0;i<NVERT;++i)
      {
	int otherId;
	const Vertex *ref=this->getVertex(i);
	if (i!=index)
	  for (otherId=0;nei->getVertex(otherId)!=ref;++otherId){}
	else otherId=oppIndex;
	order[i]=otherId;	
      }

    // Count the number of permutation to sort the indices
    int nPerm=0;
    for (int i=0;i<NVERT-1;++i)
      for (int j=i+1;j<NVERT;++j)
	{
	  if (order[i]<order[j]) 
	    {
	      nperm++;
	      std::swap(order[i],order[j]);
	    }
	}

    return ((nPerm&1)<<1)-1;
  }
  */

  int getNeighborsCount()
  {
    int count=0;
    for (int i=0;i<NNEI;i++) if (Base::neighbors[i]!=NULL) ++count;
    return count;
  } 

  LocalIndex getLocalIndex() const
  {
    return Base::localIndex;
  }
  
  GlobalIdentity getGlobalIdentity(int rank) const
  {  
    if (Base::isShadow())
      return getAsShadow()->getGlobalIdentity();
    if (Base::isGhost())
      return getAsGhost()->getGlobalIdentity();
    
    return GlobalIdentity(rank,getLocalIndex());
    //return GlobalIdentity(mpiComWorld->rank(),getLocalIndex());
  }

  GlobalIdentity getGeneration() const
  {
    return Base::generation;
  }

  int getLevel() const
  {
    return Base::generation.rank();
  }

  template <class G, int filter = predicate::filterType::Adaptive>
  bool pointIsInside(const Coord * coords, const G * geometry, 
		     int coordsAreConsistent=-1) const
  {
    return predicate::PointInSimplexT<NDIM,filter>::template test<MyType,G>
      (this,coords,geometry,coordsAreConsistent); 
  }
  
  template <class L>
  void print() const
  {
    char tmp[30];
    GlobalIdentity gid(getGlobalIdentity(GlobalIdentity::MAX_RANK));
    sprintf(tmp,"(%ld,%ld)",(long)gid.rank(),(long)gid.id());
    
    print<L>(tmp);
  }
  
  template <class L>
  void print(const std::string &s) const
  {   
    char tmp[30];
    GlobalIdentity gid(getGlobalIdentity(GlobalIdentity::MAX_RANK));
    sprintf(tmp,"(%ld,%ld)",(long)gid.rank(),(long)gid.id());
    std::string t(tmp);
    Base::template print<L>(s+t);   
  }

  // friend struct cmpPtrLocalIndexMore;
  // friend struct cmpPtrLocalIndexLess;
  struct cmpPtrLocalIndexLess
  {
    bool operator()(const MyType* a,const MyType* b)
    {
      return a->getLocalIndex()<b->getLocalIndex();
    }
  };
  struct cmpPtrLocalIndexMore
  {
    bool operator()(const MyType* a,const MyType* b)
    {
      return a->getLocalIndex()>b->getLocalIndex();
    }
  };


  // bool operator<  (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)<
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}
  // bool operator<= (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)<=
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}
  // bool operator>  (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)>
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}
  // bool operator>= (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)>=
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}
  // bool operator== (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)==
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}
  // bool operator!= (const MyType &other) const 
  // {return getGlobalIdentity(GlobalIdentity::MAX_RANK)!=
  //     other.getGlobalIdentity(GlobalIdentity::MAX_RANK);}

  /*
  bool operator<  (const MyType &other) const 
  {return Base::localIndex<other.Base::localIndex;}
  bool operator<= (const MyType &other) const 
  {return Base::localIndex<=other.Base::localIndex;}
  bool operator>  (const MyType &other) const 
  {return Base::localIndex>other.Base::localIndex;}
  bool operator>= (const MyType &other) const 
  {return Base::localIndex>=other.Base::localIndex;}
  bool operator== (const MyType &other) const 
  {return Base::localIndex==other.Base::localIndex;}
  bool operator!= (const MyType &other) const 
  {return Base::localIndex!=other.Base::localIndex;}
  */

  //protected:
public: 
  typedef typename LeafBase::Node Node;

  // splits the simplex by breaking segment seg at newVertex
  // this function only takes care of the current simplex and the caller
  // must ensure consistency with the rest of the tesselation.
  // The local ID of newSimplex should be set correctly BEFORE calling this ...
  void splitSegment(const SegmentHandle &seg, Node *newNode, Vertex *newVertex, 
		    Simplex *newSimplex, Segment **newSegment=NULL)
  {    
    LocalIndex id=newSimplex->getLocalIndex();
   
#pragma omp critical
    (*newSimplex)=(*this);
    newSimplex->setLocalIndex(id);

    /*
    if (newNode==NULL)
      {
	if (!Base::isShared())
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Second parameter 'newNode' may only be NULL for shadow or ghosts simplices ...\n");
	    exit(-1);
	  }
      }
    else LeafBase::split(newNode,newSimplex,true);
    */
 
    this->shrink(seg,newVertex,newSegment,true,newSimplex);
    newSimplex->shrink(seg,newVertex,newSegment,false,this);

    if (newNode!=NULL)
      {	
	// Be carefull : order must be the same as in split
	// => first 'this' and then 'newSimplex'
	int i1=this->getVertexIndex(newVertex);
	int i2=newSimplex->getVertexIndex(newVertex);
#pragma omp critical
	{
	  LeafBase::split(newNode,newSimplex,true);
	  LeafBase::setSplitIndices(i1,i2);
	}
      }

    // After splitting, we want to be able for every split simplex to retrieve its
    // partner, the newly created vertex and the reference simplex responsible for the
    // splitting itself. But we cannot store all that in the simplex cache so we have
    // to be a little bit sneaky ...

    // chop the pointer to reference simplex in 2 and store it in the second unsigned 
    // integer part of the cache
    this->cache.ptr=seg->getSimplex();
    newSimplex->cache.ui[1]=this->cache.ui[0];
    
    // the first and second 'char' values of the cache store the index of the partner
    // among the neighbors and the index of the new vertex respectively
    this->cache.c[0]=this->getNeighborIndex(newSimplex);
    this->cache.c[1]=this->getVertexIndex(newVertex);    
    newSimplex->cache.c[0]=newSimplex->getNeighborIndex(this);
    newSimplex->cache.c[1]=newSimplex->getVertexIndex(newVertex);

    // only the non new simplex is tagged, so that we can identify it
    this->setTaggedF();  
    newSimplex->setTaggedF(false);
    
    int newGeneration=this->getGeneration().rank()+1;    
    newSimplex->setGeneration(newGeneration,this->getGeneration().id());
    this->setGeneration(newGeneration,this->getGeneration().id());
  }

  // Retrieve the neighbors of a simplex after it was split and set the
  // boundary flags. Each neighbor is either the former neighbor before 
  // splitting or its partner ...
  // Note that the cache.c[0] value of all neighboring split simplex should be set to the 
  // index of the partner simplex in its neighbors list before calling this function
  bool fixNeighborsAfterSplitting(bool skipShared=true, bool updateVertexFlags=true)
  {
    Vertex *curV[NVERT];
    Vertex *neiV[NVERT];	
    bool ret=false;
    bool bndry=false;
    Base::getVertices(curV);
    
    for (int i=0;i<NNEI;i++)
      {
	if (Base::neighbors[i]==NULL) 
	  {
	    // shadows do not maintain boundary informations
	    // as their neighbors do not have to be set
	    // But for other simplices, that means we are on the boundary !
	    if (!Base::isShadow())
	      {
		bndry=true;
		if (updateVertexFlags)
		  {
		    for (int j=0;j<NVERT;++j)
		      if (j!=i) Base::vertices[j]->setBoundaryF();
		  }
	      }
	    continue;
	  }

	if ((skipShared)&&(Base::neighbors[i]->isShadowOrGhost()))
	  continue;


	Base::neighbors[i]->getVertices(neiV);
	int nMissed=0;
	for (int k=0;k<NVERT;k++)
	  {
	    int j=0;
	    for (;j<NVERT;j++)
	      {
		if (curV[j]==neiV[k]) break;
	      }
	    if (j==NVERT) nMissed++;
	  }
	
	
	if (nMissed!=1)
	  {
	    // If a neighbor is wrong, then it's partner should be fine ...
	    Base::neighbors[i] = Base::neighbors[i]->
	      getNeighbor(Base::neighbors[i]->cache.c[0]);
	    //Base::neighbors[i] = Base::neighbors[i]->getPartner();
	    ret =true;	    
	  }
      }
    /*
    if (nCall==0)
      {
	if (fixNeighborsAfterSplitting(1))
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Could not fix neighbors after splitting.\n");
	    exit(-1);
	  }
	else glb::console->print<LOG_DEBUG>("Fixed\n");
      }
    else
      {
	print<LOG_DEBUG>("Current: ");
	for (int i=0;i<NNEI;i++)
	  Base::neighbors[i]->template print<LOG_DEBUG>("  =>");
      }
    */

    Base::setSafeF(true);
    Base::setBoundaryF(bndry);
    
    return ret;
  }

   // retrieve the partner simplex (i.e. the other part of the split parent simplex)
  Simplex *getPartner()
  {
    // shared simplices do not belong to a tree structure
    if (Base::isShared())
      {
	// Temporary partners are stored in LeafBase::parent and are 
	// valid only until one of the partners is split again
	//return static_cast<Simplex*>(LeafBase::parent);
	return static_cast<Simplex*>(NULL);
      }
    else return static_cast<Simplex*>(LeafBase::getPartner());      
  }

  Simplex *getPartnerNotShared() const
  {
    return static_cast<Simplex*>(LeafBase::getPartner());      
  }

protected:

  template <class SU,class GSU, class SSU>
  void updateSimplicesPointers(SU &su, GSU &gsu, SSU &ssu)
  {
    for (int i=0;i<NNEI;++i)
      {
	Simplex *nei = Base::neighbors[i];
	if (nei!=NULL)
	  {
	    if (nei->isLocal())
	      Base::neighbors[i]=su(nei);
	    else if (nei->isGhost())
	      Base::neighbors[i]=gsu(static_cast<Ghost*>(nei));
	    else
	      Base::neighbors[i]=ssu(static_cast<Shadow*>(nei));
	  }
      }

    if (Base::isLocal()) LeafBase::updateParentPointer(su);
      
  }

  template <class VU,class GVU, class SVU>
  void updateVerticesPointers(VU &vu, GVU &gvu, SVU &svu)
  {
    for (int i=0;i<NVERT;++i)
      {
	Vertex *v=Base::vertices[i];
	if (v->isLocal())
	  Base::vertices[i]=vu(v);
	else if (v->isGhost())
	  Base::vertices[i]=gvu(static_cast<GhostVertex*>(v));
	else
	  Base::vertices[i]=svu(static_cast<ShadowVertex*>(v));
      }    
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
    typedef typename VU::DataPtr VP;
    typedef typename GVU::DataPtr GVP;
    typedef typename SVU::DataPtr SVP;

    typedef typename SU::DataPtr SP;
    typedef typename GSU::DataPtr GSP;
    typedef typename SSU::DataPtr SSP;

    for (int i=0;i<NVERT;++i)
      {		
	if (Base::vertices[i]!=NULL)
	  {
	    Vertex *v=vertexUpdate(Base::vertices[i]);
	    if (v==NULL) v=ghostVertexUpdate(static_cast<GVP>(Base::vertices[i]));
	    if (v==NULL) v=shadowVertexUpdate(static_cast<SVP>(Base::vertices[i]));
	    if (v==NULL)
	      {
		PRINT_SRC_INFO(LOG_ERROR);
		glb::console->print<LOG_ERROR>("when unserializing : could not update vertex pointer !\n");
		exit(-1);
	      }
	    Base::vertices[i]=v;
	  }
      }
    for (int i=0;i<NNEI;++i)
      {
	if (Base::neighbors[i]!=NULL)
	  {
	    Simplex *s=simplexUpdate(Base::neighbors[i]);
	    if (s==NULL) s=ghostSimplexUpdate(static_cast<GSP>(Base::neighbors[i]));
	    if (s==NULL) s=shadowSimplexUpdate(static_cast<SSP>(Base::neighbors[i]));
	    if (s==NULL)
	      {
		PRINT_SRC_INFO(LOG_ERROR);
		glb::console->print<LOG_ERROR>("when unserializing : could not update neighbor simplex pointer !\n");
		exit(-1);
	      }
	    Base::neighbors[i]=s;
	  }
      }
  }

  // current and neighbors simplices must be set (e.g. with setSegments() or set()) before 
  // calling this function
  // this function also sets the boundary/shared flags correctly for the simplex and its segments/vertices
  // when one ore more neigbors is NULL or is a shadow
  template <typename inputIteratorN>
  void setNeighbors(inputIteratorN &nei)
  {   
    std::copy(nei,nei+NNEI,Base::neighbors); 

    //flags |= SIMPLEX_FLAG_UNSAFE;
    Base::setSafeF(false);

    if (!Base::isSet())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Elements must be set before neighbors (e.g. with set()).\n");
	exit(-1);
      }
    
    sortNeighbors();    
    
    // Check if we are on a local / global boundary and tag accordingly
    for (int i=0;i<NNEI;i++)
      {
	if (Base::neighbors[i]==NULL) // global boundary
	  {
	    Base::setBoundaryF(true);
	    Element *el[NOPPEL];	 
	    Base::getOppositeElements(el,i); 
	    for (int j=0;j<NOPPEL;j++) el[j]->setBoundaryF(true);	    
	  }
	/*
	else if (Base::neighbors[i]->isShadow()) // local boundary
	  {
	    Base::setSharedF(true);
	    Element *el[NOPPEL];	 
	    Base::getOppositeElements(el,i); 
	    for (int j=0;j<NOPPEL;j++) el[j]->setSharedF(true);
	  }
	*/
      }

    Base::setSafeF();
  }

  // Ensure that neighbors[i] corresponds to the neighbor opposite to vertex(i)
  bool sortNeighbors() 
  {
    if (Base::isSafe()) return false;
    
    Simplex *sortedNei[NNEI];

    for (int i=0;i<NNEI;i++) sortedNei[i]=NULL;

    for (int i=0;i<NNEI;i++)
      {
	if (Base::neighbors[i]==NULL) continue;

	int id = getOppositeVertexIndex(Base::neighbors[i]);

	if (id<0)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("Simplices are not neighbors.\n");	    
	    exit(-1);
	  }
	
	sortedNei[id]=Base::neighbors[i];
      }

    std::copy(sortedNei,sortedNei+NNEI,Base::neighbors);
    
    return true;
  } 

  void setLocalIndex(LocalIndex id)
  {
    Base::localIndex = id;
  }
  
  void setGeneration(long gen,unsigned long index)
  {    
    Base::generation = GlobalIdentity(gen,index);
  }

  void setGeneration(GlobalIdentity gen)
  {    
    Base::generation = gen;
  }

  
  template <class LOG>
  bool checkConsistency(bool report=true)
  {
    bool success=true;
    int localVertexIndex=-1;
    for (int i=0;i<NVERT;++i)
      {
	Vertex *v=this->getVertex(i);
	if (v==NULL)
	  {
	    if (report)
	      glb::console->print<LOG>("NULL vertex @%d.\n",i);
	    success=false;
	  }
	else
	  {
	    if ( (v->isShadow()) && (!Base::isShadow()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Shadow vertex in non shadow simplex @%d.\n",i);
		    v->template print<LOG>();
		  }
		success=false;
	      }
	    if ( (v->isShadowOrGhost()) && (!Base::isShadowOrGhost()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Ghost/shadow vertex in local simplex @%d.\n",i);
		    v->template print<LOG>();
		  }
		success=false;
	      } 
	    if ( (!v->isShadowOrGhost()) && (Base::isShadow()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Local vertex in shadow simplex @%d.\n",i);
		    v->template print<LOG>();
		  }
		success=false;
	      } 
	    if ( (!v->isShared()) && (Base::isShadowOrGhost()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Non shared vertex in shadow/ghost simplex @%d.\n",i);
		    v->template print<LOG>();
		  }
		success=false;
	      } 
	    if (!v->isShadowOrGhost()) localVertexIndex=i;
	  }
      }

    if ((Base::isShadow())&&(localVertexIndex>=0))
      {
	if (report)
	  {
	    glb::console->print<LOG>("Shadow simplex has at least one local vertex @%d.\n",
				     localVertexIndex);
	    this->getVertex(localVertexIndex)->template print<LOG>();
	  }
	success=false;
      }

    if (!success) return success;

    for (int i=0;i<NNEI;++i)
      {
	Simplex *nei=this->getNeighbor(i);
	// FIXME: check for NULL neighbor in non boundary simplex ?
	if (nei==NULL)
	  {
	    
	  }
	else
	  {
	    if ( (Base::isShadow()) && (!nei->isShadowOrGhost()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Local neighbor in shadow simplex @%d.\n",i);
		    nei->template print<LOG>("NEI:");
		  }
		success=false;
	      }
	    if ( (!Base::isShadowOrGhost()) && (nei->isShadow()) )
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Shadow neighbor in local simplex @%d.\n",i);
		    nei->template print<LOG>("NEI:");
		  }
		success=false;
	      }
	    
	    if (nei->getNeighborIndex(this)<0)
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Neighborhood relation is not reciprocal ! @%d.\n",i);
		    nei->template print<LOG>("NEI:");
		    this->getVertex(i)->template print<LOG>();
		  }
		success=false;
	      }
	    
	    int nVert=0;
	    for (int j=0;j<NVERT;++j)
	      if (nei->getVertexIndex(this->getVertex(j))>=0) 
		nVert++;
	    if (nVert != NVERT-1)
	      {
		if (report)
		  {
		    glb::console->print<LOG>("Simplex shares %d vertices with neighbor @%d (should be %d).\n",nVert,i,NVERT-1);
		    nei->template print<LOG>("NEI:");
		  }
		success=false;
	      }
	  }
      }

    return success;
  }

private:
  // These function are used only when building the simplices, while the 
  // neighbors are not correctly set (i.e. the simplex is unsafe)

  /*
  // suppose that nei is indeed a neighbor ....
  Vertex *getOppositeVertex(Simplex *nei) const
  {       
    
    // if (Base::isSafe())
    //   {
    // 	for (int i=0;i<NNEI;i++)
    // 	  if (Base::neighbors[i]==nei) return Base::getVertex(i);
    // 	return NULL;
    //   }
    
    if (nei->isShadow()) // in that case, 'nei' knows better than me ...
      {
	return nei->getAsShadow()->getNeighborOppositeVertex(this);
      }
    
    if (!nei->isSet())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Requested neighbor is not set (call nei->set() or nei->setSegments() before this !!!).\n");
	exit(-1);
      }

    hlp::FindFirstDifference<Vertex*,NNEI> firstDiff;
    Vertex *curV[NVERT];
    Vertex *neiV[NVERT];	

    Base::getVertices(curV);
    nei->getVertices(neiV);
    int id=firstDiff(curV,neiV);

    if (id<0) return NULL; // nei is identical to the current simplex

    // check wether nei is a neighbor 
    // <=> only curV[id] is not a vertex of nei
    for (int i=id+1;i<NVERT;i++)
      {
	int j=0;
	for (;j<NVERT;j++)
	  {
	    if (curV[i]==neiV[j]) break;
	  }
	if (j==NVERT) return NULL; // not neighbors !!!
      }

    return curV[id]; // OK found  
  }
  */
  
  int getOppositeVertexIndex(Simplex *nei) const
  { 
    /*
    if (Base::isSafe())
      {
	for (int i=0;i<NNEI;i++)
	  if (Base::neighbors[i]==nei) return i;
	return -1;
      }
    */

    if (!nei->isSet())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Requested neighbor is not set (call nei->set() or nei->setSegments() before this !!!).\n");
	exit(-1);
      }

    //hlp::FindFirstDifference<Vertex*,NNEI> firstDiff;
    Vertex *curV[NVERT];
    Vertex *neiV[NVERT];	

    Base::getVertices(curV);
    nei->getVertices(neiV);
    //if ((curV[0]==curV[1])||(curV[1]==curV[2])||(curV[0]==curV[2])) exit(-1);
    int id=hlp::FindFirstDifference<Vertex*,NNEI>::find(curV,neiV);

    if (id<0) return -1; // nei is identical to the current simplex

    // check wether nei is actually a neighbor 
    // <=> only curV[id] is not a vertex of nei
    for (int i=id+1;i<NVERT;i++)
      {
	int j=0;
	for (;j<NVERT;j++)
	  {
	    if (curV[i]==neiV[j]) break;
	  }
	if (j==NVERT) return -1; // not neighbors !!!
      }

    return id; // OK found
  } 

  void swap(Ghost *other)
  {
    // only shadow and ghosts can be swapped as they do not need the tree information
    // which is NOT swappable

    Base::swap_internal(other);
    if (Base::isShadow())
      other->setGlobalIdentity(static_cast<Shadow*>(this)->getGlobalIdentity());
    else if (Base::isGhost())
      other->setGlobalIdentity(static_cast<Ghost*>(this)->getGlobalIdentity());
    else
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Only Ghost and/or Shadow simplices may be swapped\n");
	exit(-1);
      }
  }

  void swap(Shadow *other)
  {
    // only shadow and ghosts can be swapped as they do not need the tree information
    // which is NOT swappable

    Base::swap_internal(other);
    if (Base::isShadow())
      other->setGlobalIdentity(static_cast<Shadow*>(this)->getGlobalIdentity());
    else if (Base::isGhost())
      other->setGlobalIdentity(static_cast<Ghost*>(this)->getGlobalIdentity());
    else
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Only Ghost and/or Shadow simplices may be swapped\n");
	exit(-1);
      }
  }

};

/** \}*/
#include "../internal/namespace.footer"
#endif
