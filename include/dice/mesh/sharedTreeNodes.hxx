#ifndef __SHARED_TREE_NODES_HXX__
#define __SHARED_TREE_NODES_HXX__

/**
 * @file 
 * @brief  A class used to define binary tree nodes and leaves.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

template <class E, class T> class SharedTreeT;

/**
 * \namespace sharedTreeNodes 
 * \brief classes used to define binary tree nodes and leaves (see SharedTreeT).
 *
 * This version emulates virtual functions instead of actually using virtual 
 * inheritence.
 * Because there is no Vtable needed, we can store 8 more bytes per object 
 * to identify their type AND store some more info ! In particular, we store for each
 * leaf whether its got an available partner or not (i.e. the other simplex resulting 
 * from a split).  In this version, the "virtual" functions are defined in AnyT<E,T>, 
 * depending on the type as stored by the GlobalIdentity.type() value 
 * (see isLeaf, isNode and isRoot).
 *
 * \todo Using CRTP may be more elegant ...
 */
namespace sharedTreeNodes {
  
  template <class E, class T> class AnyT;
  template <class E, class T> class NodeT;
  template <class E, class T> class LeafT;
  template <class E, class T> class RootT;

  /**
   * \class AnyT 
   * \brief Base for tree nodes. Represents any type of nodes.
   */
  template <class E, class T>
  class AnyT
  {
  public:  
    typedef E Element;
    typedef T Traits;

    typedef AnyT<E,T>  Any;
    typedef NodeT<E,T> Node;
    typedef RootT<E,T> Root;
    typedef LeafT<E,T> Leaf;

    friend class SharedTreeT<E,T>;
    //friend class AnyT<E,T>;
    friend class NodeT<E,T>;
    friend class RootT<E,T>;
    friend class LeafT<E,T>;

    typedef AnyT<E,T> MyType;
    typedef AnyT<E,T> AnyNodeBase;
    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename GlobalIdentity::Value GlobalIdentityValue;
    typedef typename GlobalIdentity::Type GlobalIdentityType;
    typedef typename GlobalIdentity::Id GlobalIdentityId;

    static const GlobalIdentityType leafType = GlobalIdentity::MAX_TYPE;
    static const GlobalIdentityType nodeType = GlobalIdentity::MAX_TYPE-1;

    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {      
      writer->write(&me->parent);
      writer->write(&me->gid);
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      reader->read(&me->parent);
      reader->read(&me->gid);
    }

    bool isRoot() const 
    {
      return gid.type() < nodeType;
    }
  
    bool isNode() const 
    {
      return gid.type() == nodeType;
    }

    bool isLeaf() const 
    {
      return gid.type() == leafType;
    }
    
    MyType *getOther(const MyType *cur) 
    {
      MyType *o=NULL;
      if (isNode())
	{
	  MyType** children=static_cast<Node*>(this)->getChildren();
	  if (children[0]==cur) 
	    o=children[1];
	  else if (children[1]==cur) 
	    o=children[0];
	}
      return o;
    }

    Root *getRoot()
    {
      if (isRoot())
	return static_cast<Root*>(this);
      else
	return parent->getRoot();
    }

    const Root *getConstRoot() const
    {
      if (isRoot())
	return static_cast<const Root*>(this);
      else
	return parent->getRoot();
    }

    int getWeight() const
    {
      if (isRoot())
	return static_cast<const Root*>(this)->weight;
      else
	return parent->getWeight();
    }

    int getDepth() const 
    {
      if (isRoot())
	return -1;
      else
	return 1+parent->getDepth();
    }

    MyType *getParent()
    {
      if (!isRoot())
	return parent;
      else
	return NULL;
    }    
    /*
    MyType *getChild(int which)
    {
      if (isRoot())
	return static_cast<Root*>(this)->getChild();
      else if (isNode())
	return static_cast<Node*>(this)->getChild(which);
      else return NULL;	
    }
    */
    int addChild(MyType *c)
    {
      if (isRoot())
	{
	  if (static_cast<Root*>(this)->getChild() == NULL)
	    {
	      static_cast<Root*>(this)->setChild(c);
	      c->setParent(this);
	      return 0;
	    }
	}
      else if (isNode())
	{
	  MyType** children = static_cast<Node*>(this)->getChildren();
	  if (children[0]==NULL)
	    {
	      children[0]=c;
	      c->setParent(this);
	      return 0;
	    }
	  else if (children[1]==NULL)
	    {
	      children[1]=c;
	      if (children[0]->isLeaf()&&children[1]->isLeaf())
		{
		  static_cast<Leaf*>(children[0])->
		    setSplitIndex(static_cast<const Node*>(this)->getChildIndex(0));
		  static_cast<Leaf*>(children[1])->
		    setSplitIndex(static_cast<const Node*>(this)->getChildIndex(1));
		  
		  //static_cast<Leaf*>(children[0])->haveDirectPartner(true);
		  //static_cast<Leaf*>(children[1])->haveDirectPartner(true);
		}	
	      c->setParent(this);
	      return 1;
	    }	  
	}
      
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("Cannot add any more child (full)");
      exit(-1);

      return -1;
    }

  protected:
    MyType *parent;
    GlobalIdentity gid;

    void setParent(MyType *p)
    {
      if (!isRoot()) parent=p;
      else
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Should never be called from a ROOT node.\n");
	  exit(-1);
	}
    }
    
    template <class IT>
    long buildCompactTree(std::vector<IT> &compactTree, char &index, 
			  std::vector<IT> &nodeData)
    {
      if (isLeaf())
	{
	  if (index==0) compactTree.push_back(0);   
	  if ((++index) == sizeof(IT)*8) index=0;   
     
	  return 1;
	}
      else if (isNode())
	{
	  if (index==0) compactTree.push_back(0);      
    
	  compactTree.back() |= (1ul << index);
	  if ((++index) == sizeof(IT)*8) index=0;
	  nodeData.push_back((IT)gid.id());
	  Node* me=static_cast<Node*>(this);
	  long n=me->children[0]->buildCompactTree(compactTree,index,nodeData);
	  n+=me->children[1]->buildCompactTree(compactTree,index,nodeData);
	  
	  return n;
	}
      else
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Should never be called from a ROOT node.\n");
	  exit(-1);	  
	}  
      return 0;
    }
    
    template <class MP>
    Node* getLeavesRecycleNodes(std::vector<Element*> &leaves, MP &nodePool)
    {
      if (isLeaf())
	{
	  leaves.push_back(static_cast<E*>(this));
	  return NULL;
	}
      else if (isNode())
	{
	  Node* res;  
	  MyType** children=static_cast<Node*>(this)->getChildren();
	  res = children[0]->getLeavesRecycleNodes(leaves,nodePool);
	  if (res!=NULL) nodePool.recycle(res);

	  res = children[1]->getLeavesRecycleNodes(leaves,nodePool);
	  if (res!=NULL) nodePool.recycle(res);

	  return static_cast<Node*>(this);
	}
      else
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Should never be called from a ROOT node.\n");
	  exit(-1);
	}
      return NULL;
    }

  protected:
    // This class may not be constructed ! It's meant to be inherited ...
    AnyT() {}
  };

  /**
   * \class LeafT 
   * \brief A Leaf type node.
   * Meaning of gid.id():
   * 0 -> not partner
   * 1 -> partner present, don't know the index
   * 2+ -> partner present, index = gid.id()-2
   */
  
  template <class E, class T>
  class LeafT : public AnyT<E,T>
  {
  public:  
    typedef E Element;
    typedef T Traits;

    typedef NodeT<E,T> Node;
    typedef RootT<E,T> Root;
    typedef LeafT<E,T> Leaf;

    friend class SharedTreeT<E,T>;
    friend class NodeT<E,T>;
    friend class RootT<E,T>;
    //friend class LeafT<E,T>;
    friend class AnyT<E,T>;

    typedef LeafT<E,T> MyType;
    typedef AnyT<E,T>  Base;
    typedef AnyT<E,T>  AnyNodeBase;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename GlobalIdentity::Value GlobalIdentityValue;

    LeafT()
    {
      Base::parent=NULL;
      Base::gid=GlobalIdentity(Base::leafType,0);
    }    

    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {      
      Base::selfSerialize(static_cast<const Base*>(me), writer);
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      Base::selfUnSerialize(static_cast<Base*>(me), reader);
    }

    bool haveDirectPartner() const 
    {
      return (bool)Base::gid.id();
    }

    int getSplitIndex() const
    {
      return ((int)Base::gid.id())-2;

      // if (Base::parent->isNode())
      // 	return static_cast<Node*>(Base::parent)->getChildIndex(this);
	
      // return -1;
    }
    
    template <class SU>    
    void updateParentPointer(SU &su)
    {
      if (Base::parent->isRoot())
	static_cast<Root*>(Base::parent)->setChild(this);
      else
	{
	  Node *p=static_cast<Node*>(Base::parent);
	  if (p->getChild(0)->isLeaf())
	    {
	      LeafT *l0=static_cast<LeafT*>(p->getChild(0));
	      Element *e0=static_cast<Element*>(l0);
	      l0=static_cast<LeafT*>(su(e0));

	      if ( l0 == this )
		p->replaceChild(p->getChild(0),this);
	      else
		p->replaceChild(p->getChild(1),this);
	    }
	  else p->replaceChild(p->getChild(1),this);
	  
	}
      //setChild(c);
    }
   
  protected:

    void haveDirectPartner(bool val)
    {
      if (val)
	{
	  // if we laready have an index, preserve it
	  if (Base::gid.id()<1)
	    Base::gid.setId(1);
	}
      else Base::gid.setId(0);	
	
    }

    void setSplitIndex(int index)
    {
      Base::gid.setId(2+index);
    }
    
    MyType* getPartner() const
    {
      if (haveDirectPartner())
	return static_cast<MyType*>(Base::parent->getOther(this));
      else
	return NULL;
    }

    // void setParent(Base *b)
    // {
    //   Base::parent=b;
    // }
    
    
    void split(Node *node, Leaf *leaf, bool preventWeightIncrease=false)
    {
      if (node==NULL) return;
      
      MyType *p=getPartner();
      if (p!=NULL)
	p->haveDirectPartner(false);

      if (Base::parent->isRoot())
	static_cast<Root*>(Base::parent)->setChild(node);
      else if (Base::parent->isNode())
	static_cast<Node*>(Base::parent)->replaceChild(this,node);	  
	         
      node->setParent(Base::parent);
      node->setChildren(this,leaf);  
   
      leaf->setParent(node);  
      leaf->haveDirectPartner(true);  
      
      Base::setParent(node);
      haveDirectPartner(true);

      if (!preventWeightIncrease)
	{
	  //Base::getRoot()->increaseWeight();
	  leaf->getRoot()->increaseWeight();
	}
    }
    
    void setSplitIndices(unsigned char i1, unsigned char i2)
    {
      static_cast<Node*>(Base::parent)->setChildrenIndex(i1,i2);  
    }
   
    // one of the leaves is returned for removal
    template <class NP>
    Element *merge(NP &nodePool)
    {
      MyType *partner=getPartner();
      if (partner==NULL)
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Leaf has no partner, cannot merge it!\n");
	  exit(-1);
	} 
      
      Node *parent = static_cast<Node*>(Base::parent);
      Base *grandParent = parent->getParent();      

      Leaf *keep = static_cast<Leaf*>(parent->getChild(0));
      Leaf *rm = static_cast<Leaf*>(parent->getChild(1));

      if (grandParent->isRoot())
	{
	  static_cast<Root*>(grandParent)->setChild(keep);
	  static_cast<Root*>(grandParent)->decreaseWeight();
	  keep->setParent(grandParent);
	  keep->haveDirectPartner(false);
	}
      else
	{	  
	  static_cast<Node*>(grandParent)->replaceChild(parent,keep);
	  static_cast<Node*>(grandParent)->getRoot()->decreaseWeight();
	  keep->setParent(grandParent);
	  Base *other = grandParent->getOther(keep);
	  if (other->isLeaf())
	    {	 
	      keep->
		setSplitIndex(static_cast<Node*>(grandParent)->getChildIndex(keep));
	      static_cast<Leaf*>(other)->
		setSplitIndex(static_cast<Node*>(grandParent)->getChildIndex(other));
	    }
	  else keep->haveDirectPartner(false);
	}
            
      nodePool.recycle(parent);
      return static_cast<Element*>(rm);
    }
  };

  /**
   * \class NodeT 
   * \brief Any node that is not a root or a leaf.
   */
  template <class E, class T>
  class NodeT : public AnyT<E, T>
  {
  public: 
    typedef E Element;
    typedef T Traits;

    typedef NodeT<E,T> Node;
    typedef RootT<E,T> Root;
    typedef LeafT<E,T> Leaf;

    friend class SharedTreeT<E,T>;
    //friend class NodeT<E,T>;
    friend class RootT<E,T>;
    friend class LeafT<E,T>;
    friend class AnyT<E,T>;    

    typedef NodeT<E,T> MyType;
    typedef AnyT<E,T>  Base;
    typedef AnyT<E,T>  AnyNodeBase;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename GlobalIdentity::Value GlobalIdentityValue;

    NodeT()
    {
      Base::parent=NULL;
      Base::gid = GlobalIdentity(Base::nodeType,0);
      children[0]=NULL;
      children[1]=NULL;
    }

    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {     
      writer->write(me->children,2);
      Base::selfSerialize(static_cast<const Base*>(me), writer);
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      reader->read(me->children,2);
      Base::selfUnSerialize(static_cast<Base*>(me), reader);
    }
    
    Base* getChild(int i)
    {
      return children[i];
    }

    Base** getChildren()
    {
      return children;
    }
   
  protected:
    Base* children[2];

    void setChildrenIndex(unsigned char i1, unsigned char i2)
    {
      Base::gid.setId((static_cast<long>(i1)<<8)+
		      static_cast<long>(i2));
      if ((children[0]->isLeaf())&&(children[1]->isLeaf()))
	{
	  static_cast<Leaf*>(children[0])->setSplitIndex(i1);
	  static_cast<Leaf*>(children[1])->setSplitIndex(i2);
	}
    }

    void setChildrenIndex(long val)
    {
      Base::gid.setId(val);
    }

    int getChildIndex(int which) const
    {
      if (which==0)
	return Base::gid.id()>>8;
      else
	return Base::gid.id()&((1<<8)-1);
    }

    int getChildIndex(Base *which) const
    {
      if (which==children[0])
	return Base::gid.id()>>8;
      else if (which==children[1])
	return Base::gid.id()&((1<<8)-1);
      return -1;
    }

    void replaceChild(Base *curChild,Base *newChild)
    {
      if (children[0] == curChild)
	children[0]=newChild;
      else if (children[1] == curChild)
	children[1]=newChild;
    }

    void setChildren(Base *a,Base *b)
    {
      children[0]=a;
      children[1]=b;
    }

    template <class NU, class EU>
    void updateAfterUnserialized(const NU  &nodeUpdater,
				 const EU  &elementUpdater,
				 bool swap)
    {
      for (int i=0;i<2;++i)
	{
	  Base *c = children[i];
	  if (c!=NULL)
	    {
	      c=nodeUpdater(children[i]);
	      if (c==NULL) // leaf
		{
		  c=elementUpdater(children[i]);
		  if (c==NULL) 
		    {
		      PRINT_SRC_INFO(LOG_ERROR);
		      glb::console->print<LOG_ERROR>("when unserializing : could not update node's child[%d].\n",i);
		      exit(-1);
		    }
		  children[i]=c;
		  //if (!c->isLeaf()) exit(0);
		  Leaf *leaf = static_cast<Leaf*>(c);		  
		  leaf->setParent(this);
		}
	      else // node
		{		  	  
		  children[i]=c;
		  //if (!c->isNode()) exit(0);
		  Node *node = static_cast<Node*>(c);	
		  node->setParent(this);
		  node->updateAfterUnserialized(nodeUpdater,elementUpdater,swap);
		}	  
	    }
	}
    }

  };

  /**
   * \class RootT 
   * \brief Root type nodes.
   * NB: in Roots, child is stored in the Base::parent field !!!
   * This may not be intuitive, but spares 8 bytes per root ...
   */
  
  template <class E, class T>
  class RootT : public AnyT<E,T>
  {
  public:
    typedef E Element;
    typedef T Traits;
  
    typedef NodeT<E,T> Node;
    typedef RootT<E,T> Root;
    typedef LeafT<E,T> Leaf;

    friend class SharedTreeT<E,T>;
    friend class NodeT<E,T>;
    //friend class RootT<E,T>;
    friend class LeafT<E,T>;
    friend class AnyT<E,T>;

    typedef RootT<E,T> MyType;
    typedef AnyT<E,T>  Base;
    typedef AnyT<E,T>  AnyNodeBase;

    typedef typename T::GlobalIdentity GlobalIdentity;
    typedef typename GlobalIdentity::Value GlobalIdentityValue;
    typedef typename GlobalIdentity::Type GlobalIdentityType;
    typedef typename GlobalIdentity::Id GlobalIdentityId;

    static const int NNEI = T::NDIM+1;

    RootT()
    {
      Base::parent=NULL;
      Base::gid=GlobalIdentity(0,0);
      weight=0;
      std::fill(neighbors,neighbors+NNEI,static_cast<MyType *>(NULL));
    }   

    template <class W>
    static void selfSerialize(const MyType *me, W *writer)
    {     
      writer->write(me->neighbors,NNEI);
      writer->write(&me->weight);      
      Base::selfSerialize(static_cast<const Base*>(me), writer);
    }

    template <class R>
    static void selfUnSerialize(MyType *me, R *reader)
    {
      reader->read(me->neighbors,NNEI);
      reader->read(&me->weight);      
      Base::selfUnSerialize(static_cast<Base*>(me), reader);
    }

    void increaseWeight() {++weight;}
    void decreaseWeight() {--weight;}

    Base* getChild()
    {
      return Base::parent;
    }

    template <typename OutputIterator>
    int getNeighborsGlobalIdentity(OutputIterator out) const
    {
      int count=0;
      for (int i=0;i<NNEI;i++)
	{
	  if (neighbors[i]==NULL) continue;
	  (*out)=neighbors[i]->gid.get();
	  ++out;
	  ++count;
	}
      return count;
    }

    // get all the vertices except 'ref'
    template <typename OutputIterator>
    int getNeighbors(OutputIterator out) const
    {
      int count=0;
      for (int i=0;i<NNEI;i++)
	{
	  if (neighbors[i]==NULL) continue;
	  (*out)=neighbors[i];
	  ++out;
	  ++count;
	}
      return count;
    }

    GlobalIdentity getGlobalIdentity() const
    {
      return Base::gid;
    }

    GlobalIdentityId getLocalIndex() const
    {
      return Base::gid.id();
    }


    // Base* getParent()
    // {
    //   return NULL;
    // }

  protected:
    MyType *neighbors[NNEI]; // NB: NULL neighbors are stored at the end
    int weight;

    void setWeight(int w)
    {
      weight=w;
    }

    void setGlobalId(GlobalIdentity g)
    {
      Base::gid = g;
    }

    void setChild(Base *c)
    {
      Base::parent=c;
    }

    // void setParent(Base *c)
    // {
    //   PRINT_SRC_INFO(LOG_ERROR);
    //   glb::console->print<LOG_ERROR>("Should never be called from a ROOT node.\n");
    //   exit(-1);	  
    // }

    void setGlobalIdentity(GlobalIdentity g)
    {
      Base::gid=g;
    }
    /*
    template <class EU>
    void updateLeafPointers(const EU  &elementUpdater,bool swap)
    {
     
      Base* c = getChild();
      if (c!=NULL)
	{
	  if ()
	    c = nodeUpdater(getChild());
	  if (c==NULL) // Leaf
	    {
	      c = elementUpdater(getChild());
	      if (c==NULL)
		{
		  PRINT_SRC_INFO(LOG_ERROR);
		  glb::console->print<LOG_ERROR>("when unserializing : could not update root's child.\n");
		  exit(-1);
		}
	      setChild(c);
	      //if (!c->isLeaf()) exit(0);
	      Leaf *leaf = static_cast<Leaf*>(c);
	      leaf->setParent(this);
	    }
	  else // Node
	    {
	      setChild(c);
	      Node *node = static_cast<Node*>(c);
	      //if (!c->isNode()) exit(0);
	      node->setParent(this);
	      node->updateAfterUnserialized(nodeUpdater,elementUpdater,swap);
	    }
	}
    
    }
    */
    template <class RU, class SRU, class NU, class EU>
    void updateAfterUnserialized(const RU  &rootUpdater,
				 const SRU &shadowRootUpdater,
				 const NU  &nodeUpdater,
				 const EU  &elementUpdater,
				 bool recursive, bool swap)
    {
      for (int i=0;i<NNEI;++i)
	{
	  if (neighbors[i]!=NULL)
	    {
	      MyType *n = rootUpdater(neighbors[i]);
	      if (n==NULL)
		{
		  n = shadowRootUpdater(neighbors[i]);
		  if (n==NULL)
		    {
		      PRINT_SRC_INFO(LOG_ERROR);
		      glb::console->print<LOG_ERROR>("when unserializing : could not update root neighbor pointer !\n");
		      exit(-1);
		    }
		}
	      neighbors[i]=n;
	    }
	}
      
      if (recursive)
	{
	  Base* c = getChild();
	  if (c!=NULL)
	    {
	      c = nodeUpdater(getChild());
	      if (c==NULL) // Leaf
		{
		  c = elementUpdater(getChild());
		  if (c==NULL)
		    {
		      PRINT_SRC_INFO(LOG_ERROR);
		      glb::console->print<LOG_ERROR>("when unserializing : could not update root's child.\n");
		      exit(-1);
		    }
		  setChild(c);
		  //if (!c->isLeaf()) exit(0);
		  Leaf *leaf = static_cast<Leaf*>(c);
		  leaf->setParent(this);
		}
	      else // Node
		{
		  setChild(c);
		  Node *node = static_cast<Node*>(c);
		  //if (!c->isNode()) exit(0);
		  node->setParent(this);
		  node->updateAfterUnserialized(nodeUpdater,elementUpdater,swap);
		}
	    }
    	}
    }
    
  };

}

/** \}*/
#include "../internal/namespace.footer"
#endif
