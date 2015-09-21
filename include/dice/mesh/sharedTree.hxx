#ifndef __SHARED_TREE_HXX__
#define __SHARED_TREE_HXX__

#include <vector>
#include <string>

#include "../dice_globals.hxx"

#include "../mesh/mpiCellDataExchange.hxx"
#include "../mesh/mpiStructs/mpiExchangeStruct_sharedTreeRoot.hxx"

#include "../tools/helpers/find_unordered_maps.hxx"
#include "../tools/memory/memoryPool.hxx"
#include "../tools/memory/iterableMemoryPool.hxx"
#include "../tools/IO/paramsParser.hxx"

#include "../partition/parmetisParams.hxx"
#include "../partition/partitioner.hxx"
#include "../partition/partitionType.hxx"

/**
 * @file 
 * @brief  A class used to define a MPI shared binary tree used to handle cells. This class
 * is not designed to be used on its own and should be inherited
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/** 
 * \class SharedTreeT
 * \brief  A class used to define a MPI shared binary tree used to handle cells. This class
 * is not designed to be used on its own and should be inherited. In particular, this is 
 * used to handle the cells of an unstructured mesh, storing a tree of refined simplices 
 * in order to allow coarsening of the mesh .
 * Nodes of the tree are defined in sharedTreeNodes.
 */

template <class E, class T>
class SharedTreeT
{
public:
  typedef E Element;
  typedef T Traits;
  typedef SharedTreeT<E,T> MyType;

  static const int NNEI = Element::NNEI;
  static const int NDIM = Element::NDIM;

  static std::string classHeader() {return "shared_tree";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  enum Mode {UNSTRUCTURED, NETWORK}; 

  typedef typename E::LeafBase::Node Node;
  typedef typename E::LeafBase::Root Root;
  typedef typename E::LeafBase::Leaf Leaf;
   
  typedef typename Leaf::AnyNodeBase AnyNodeBase;

  typedef typename Traits::LocalIndex  LocalIndex;
  typedef typename Traits::GlobalIndex GlobalIndex;
  typedef typename Traits::GlobalIdentity GlobalIdentity;
  typedef typename GlobalIdentity::Value GlobalIdentityValue;
  
  typedef MemoryPoolT<Node> NodePool;
  typedef IterableMemoryPoolT<Root> RootPool;
  
  typedef typename RootPool::iterator network_iterator;

protected:

  SharedTreeT(MpiCommunication *com, double allocFactor=0.20):
    mpiCom(com),
    nodePool("Node",allocFactor),
    rootPool("Root",allocFactor),
    shadowRootPool("ShadowRoot",allocFactor)
  {
    curMode=NETWORK;
    mpiType_sharedTreeRoot=MpiExchg_SharedTreeRoot::MpiStruct::createMpiStructType();
    #pragma omp critical
    { 
      mpiTagsStart_repart = mpiCom->reserveTags(5);
      mpiTagsStart_general = mpiCom->reserveTags(5);
    }
  }

  /*
  template <class R>
  SharedTreeT(MpiCommunication *com, ParamsParser *parser, R *reader, 
	      double allocFactor=0.20):
    mpiCom(com),
    nodePool("Node",allocFactor),
    rootPool("Root",allocFactor),
    shadowRootPool("ShadowRoot",allocFactor)
  {      
    curMode=NETWORK;
    mpiType_sharedTreeRoot=MpiExchg_SharedTreeRoot::MpiStruct::createMpiStructType();
    #pragma omp critical
    { 
      mpiTagsStart_repart = mpiCom->reserveTags(5);
      mpiTagsStart_general = mpiCom->reserveTags(5);
    }
  }
*/
  virtual ~SharedTreeT()
  {
    //MPI_Type_free(&mpiType_sharedTreeRoot);
  }

  void setAllocFactor(double factor)
  {
    nodePool.setAllocFactor(factor);
    rootPool.setAllocFactor(factor);
    shadowRootPool.setAllocFactor(factor);
  }

  template <class W>
  void write(W *writer)
  {
    writer->writeHeader(classHeader(),classVersion()); 
 
    writer->write(&nLocalRootNodes);
    writer->write(&nRootNodesCum);
    writer->write(&loadImbalanceFactor);
    
    // write memory pools here 
    nodePool.serialize(writer);
    rootPool.serialize(writer);
    shadowRootPool.serialize(writer);
  }

  template <class EPU>
  void defrag(const EPU &epu)
  {
    typedef typename NodePool::UnserializedPointerUpdater NodePointerUpdater;
    typedef typename RootPool::UnserializedPointerUpdater RootPointerUpdater;
    
    NodePointerUpdater npu  = nodePool.defrag();
    RootPointerUpdater rpu  = rootPool.defrag();
    RootPointerUpdater srpu = shadowRootPool.defrag();

    bool swap = false;
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const network_iterator it_end=networkEnd();
	for (network_iterator it=networkBegin(i,glb::num_omp_threads);
	     it!=it_end;++it)	
	  {
	    (*it)->updateAfterUnserialized(*rpu,*srpu,*npu,*epu,true,swap);
	  }
	const network_iterator sit_end=shadowNetworkEnd();
	for (network_iterator it=shadowNetworkBegin(i,glb::num_omp_threads);
	     it!=sit_end;++it)	
	  {
	    (*it)->updateAfterUnserialized(*rpu,*srpu,*npu,*epu,false,swap);	    
	  }	
      }
  }

  template <class R, class EPU>
  void build(R *reader, const EPU &epu)
  {
    if (reader == NULL) return ;
    
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);

    glb::console->printFlush<LOG_INFO>("Unserializing binary tree ... ");

    reader->read(&nLocalRootNodes);
    reader->read(&nRootNodesCum);
    reader->read(&loadImbalanceFactor);    
    
    typedef typename NodePool::UnserializedPointerUpdater NodePointerUpdater;
    typedef typename RootPool::UnserializedPointerUpdater RootPointerUpdater;
    
    // restore the memory pools here
    NodePointerUpdater npu  = nodePool.unSerialize(reader);
    RootPointerUpdater rpu  = rootPool.unSerialize(reader);
    RootPointerUpdater srpu = shadowRootPool.unSerialize(reader);

    glb::console->printFlush<LOG_INFO>("done.\n");
    glb::console->printFlush<LOG_INFO>("Updating tree node pointers ... ");

    bool swap = reader->getNeedSwap();
    // FIXME: This is probably OK in parallel ... but it's unchecked yet !
    //#pragma omp parallel for
    //for (long i=0;i<glb::num_omp_threads;i++)
    #pragma omp parallel num_threads(glb::num_omp_threads)
      {
	int i=omp_get_thread_num();
	const network_iterator it_end=networkEnd();
	for (network_iterator it=networkBegin(i,glb::num_omp_threads);
	     it!=it_end;++it)	
	  {
	    (*it)->updateAfterUnserialized(*rpu,*srpu,*npu,*epu,true,swap);	  
	  }
	const network_iterator sit_end=shadowNetworkEnd();
	for (network_iterator it=shadowNetworkBegin(i,glb::num_omp_threads);
	     it!=sit_end;++it)	
	  {
	    (*it)->updateAfterUnserialized(*rpu,*srpu,*npu,*epu,false,swap);	    
	  }	
      }
    glb::console->printFlush<LOG_INFO>("done.\n");
  }
  
  Mode getMode() {return curMode;}

  void setMode(Mode m)
  {
    if (curMode!=m)
      {

      }
    curMode=m;
  }

  network_iterator networkBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return rootPool.begin();
    else
      return rootPool.begin(delta,stride);
  }
  
  network_iterator networkEnd()
  {
    return rootPool.end();
  }

  network_iterator shadowNetworkBegin(int delta=0, int stride=1)
  {
    if (stride<2)
      return shadowRootPool.begin();
    else
      return shadowRootPool.begin(delta,stride);
  }
  
  network_iterator shadowNetworkEnd()
  {
    return shadowRootPool.end();
  }

protected:  
  typedef mpiExchangeStruct::SharedTreeRootT<Element,Root> MpiExchg_SharedTreeRoot;
  typedef typename Partitioner::Index PartitionerIndex;  
  typedef typename my_dense_hash<GlobalIdentityValue,Root *>::type DenseHash;
  typedef typename DenseHash::iterator DenseHash_it;  
  typedef typename my_dense_hash<GlobalIdentityValue,GlobalIdentityValue>::type DenseHashGidGid;
  typedef typename DenseHashGidGid::iterator DenseHashGidGid_it;  

  Mode curMode;
  MpiCommunication *mpiCom;
  NodePool nodePool;
  RootPool rootPool;
  RootPool shadowRootPool;  
  
  std::vector<unsigned long> nLocalRootNodes;
  std::vector<unsigned long> nRootNodesCum;
  double loadImbalanceFactor;

  MpiDataType mpiType_sharedTreeRoot;
  int mpiTagsStart_repart;
  int mpiTagsStart_general;

  // can only call once !
  template <class InputIterator>
  void initRoot(InputIterator start,InputIterator stop, long N=-1)
  {
    if (N<0) N=std::distance(start,stop); 
    const int myRank=mpiCom->rank();
   
    if (curMode == NETWORK)
      {
	rootPool.reserve(N);
	//printf("max load factor : %f\n",rootHash.max_load_factor());exit(0);
	if (mpiCom->size()>1)
	  {
	    long nShadowEstimate = pow(N/mpiCom->size(),
				       ((double)NDIM-1)/((double)NDIM))*3.0;	    
	    shadowRootPool.reserve(nShadowEstimate);
	    //shadowRootHash.rehash(nShadowEstimate*1.1);
	  }

	// FIXME: set a more conservative value for max_load_factor()?
	// default is 0.5, but 0.7/0.8 is probably good enough
	//rootHash.rehash(N*(1.1/rootHash.max_load_factor()));

	InputIterator it=start;
	while(it!=stop)
	  {
	    Root *root;
	    rootPool.pop(&root);
	    root->setChild(*it);
	    root->getChild()->setParent(root);
	    root->setGlobalIdentity(it->getGlobalIdentity(myRank));
	    root->setWeight(1);

	    // root->child = static_cast<Leaf*>(*it);
	    // root->child->parent = root;
	    // root->id = it->getGlobalIdentity(myRank);	    
	    // root->weight=1;
	    
	    ++it;
	  }
	
	InputIterator it2=start;
	while(it2!=stop)
	  {
	    Element* nei[NNEI];
	    int nnei=it2->getNeighbors(nei);
	    Leaf *leaf=static_cast<Leaf*>(*it2);
	    Root *root = static_cast<Root*>(leaf->parent);

	    int j=0;
	    for (int i=0;i<nnei;i++)
	      {
		if (nei[i]==NULL) continue;
		// NB: two roots with the same remote neighbor will each have
		// a different copy of it's shadow, and that's OK ...		
		// Besides, we don't really need the tree for shadow roots,
		// but we initialize it here anyway for convenience ... 
		// It will not be correclty updated though.
		if ((nei[i]->isShadow())||(nei[i]->isGhost()))
		  {
		    Root *sharedRoot;
		    shadowRootPool.pop(&sharedRoot);
		    sharedRoot->setChild(nei[i]);
		    sharedRoot->getChild()->setParent(sharedRoot);
		    sharedRoot->setGlobalIdentity(nei[i]->getGlobalIdentity(myRank));
		    
		    // sharedRoot->child = static_cast<Leaf*>(nei[i]);
		    // sharedRoot->child->parent = sharedRoot;
		    // sharedRoot->id = nei[i]->getGlobalIdentity(myRank);
		  }

		Leaf *neiLeaf =static_cast<Leaf*>(nei[i]);
		//NB: non NULL neighbors must be stored first	
		root->neighbors[j++] = static_cast<Root*>(neiLeaf->parent); 		
	      }
	    if (j<NNEI) for (;j<NNEI;) root->neighbors[j++]=NULL;
	    ++it2;
	  }
      }
    else
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("UNSTRUCTURED mode is not implemented (and will probably never be ...)\n");
	glb::console->print<LOG_ERROR>("This function should never have been called.\n");
	exit(-1);

	rootPool.reserve(N);
	InputIterator it=start;
	while(it!=stop)
	  {
	    Root *root;
	    rootPool.pop(&root);
	    root->setChild(*it);
	    root->getChild()->setParent(root);
	    // root->child = static_cast<Leaf*>(*it);
	    // root->child->parent = root;
	    ++it;
	  }
      }   
  }

  void popRoot(Root **root)
  {
    rootPool.pop(root);
  }

  template <class InputIterator>
  void popRoots(InputIterator begin, InputIterator end)
  {
    for (InputIterator it = begin;it!=end;++it)
      rootPool.pop(&(*it));
  }

  unsigned long getNRootNodes()
  {
    if (curMode==NETWORK)
      return rootPool.getUsedCount();
    else
      return rootPool.getUsedCount();
  }

  double getLoadImbalanceFactor()
  {
    if (curMode==NETWORK)
      return loadImbalanceFactor;
    else return -1;
  }
  
  template <class ET>
  GlobalIndex getGlobalIndex(ET *element)
  {
    return globalIdentityToGlobalIndex(element->getGlobalIdentity());
  }

  GlobalIndex globalIdentityToGlobalIndex(GlobalIdentity identity)
  {
    return nRootNodesCum[identity.rank()] + identity.id();
  }
  
  double updateLoadImbalanceFactor()
  {
    unsigned long min=nLocalRootNodes[0];
    unsigned long max=nLocalRootNodes[0];
    
    for (int i=1;i<mpiCom->size();i++)
      {
	if (nLocalRootNodes[i]>max) max=nLocalRootNodes[i];
	if (nLocalRootNodes[i]<min) min=nLocalRootNodes[i];
      }
    loadImbalanceFactor = double(max)/double(min);
    return loadImbalanceFactor;
  }

  void updateRootNodesCount()
  {       
    unsigned long nNodes = getNRootNodes();
    nRootNodesCum.assign(mpiCom->size()+1,0);
    nRootNodesCum[mpiCom->rank()+1]=nNodes;	
         
    mpiCom->Allreduce_inplace(nRootNodesCum,MPI_SUM);
    nLocalRootNodes.assign(++(nRootNodesCum.begin()),nRootNodesCum.end());
    for (int i=0;i<mpiCom->size();i++)
      nRootNodesCum[i+1]+=nRootNodesCum[i];
    
    if (curMode!=NETWORK) return;
    
    updateLoadImbalanceFactor();    
  }

  std::vector<Root *> getRootNodesArray()
  {
    std::vector<Root *> result(getNRootNodes());
#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {
	const network_iterator it_end=networkEnd();
	for (network_iterator it=networkBegin(i,glb::num_omp_threads);it!=it_end;++it)	
	  {
	    result[it->getLocalIndex()]=*it;
	  }
      }
    return result;
  }
 
  bool generateParmetisGraph(ParmetisParams &p, 
			     RefinePartitionType type, 
			     double tolerance=1.05,
			     double weightPerCell=1.0)
  {
    //typedef typename ParmetisParams::Index Index;
    typedef typename ParmetisParams::Float PPFloat;

    // Weigth are integers so we need to scale them by a large enough factor,
    // but not too large so that we do not overflow an int capacity ...
    // => factor 100 means we have 2 digits precision ...
    auto minMaxSum=mpiCom->minMaxSum(weightPerCell);
    double pwFactor=1.0;
    // if (min == max) then all the weights can remain =1.0
    if (minMaxSum.first.first != minMaxSum.first.second)
      {	
	double avg=minMaxSum.second / mpiCom->size();
	pwFactor =  std::max(1.0,100.0 * (weightPerCell / avg));	
      }
    
    if (curMode != NETWORK) return false;

    updateRootNodesCount();
    p.freeData();
    p.setDefault();

    const unsigned long nRootNodes=getNRootNodes();
    p.owner=true;
    p.wgtflag=2;
    p.setTolerance(tolerance);
    
    p.vtxdist=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*nRootNodesCum.size());
    for (unsigned long i=0;i<nRootNodesCum.size();i++)
      p.vtxdist[i]=nRootNodesCum[i];   

    // PH partitionner only needs positions, not a graph ...
    if (type == RefinePartitionTypeV::PH)
      {
	p.xyz=(PPFloat*)malloc(sizeof(PPFloat)*nRootNodes*NDIM);
	p.nDims=NDIM;

	// In that case, we need positions in p.xyz
#pragma omp parallel for
	for (long i=0;i<glb::num_omp_threads;i++)
	  {	
	    const network_iterator it_end=networkEnd();
	    for (network_iterator it=networkBegin(i,glb::num_omp_threads);it!=it_end;++it)
	      {	 
		const unsigned long id = static_cast<unsigned long>(it->getLocalIndex())*NDIM;
		AnyNodeBase *any=it->getChild();

		while (!any->isLeaf()) any=static_cast<Node*>(any)->getChild(0);
		Leaf *leaf = static_cast<Leaf*>(any);
		//Element *el = static_cast<Element*>(leaf);
		//typedef typename Element::Coord Coord;
		typedef typename Element::Vertex Vertex;
		Vertex *v=static_cast<Element*>(leaf)->getVertex(0);
		for (int j=0;j<NDIM;++j) p.xyz[id+j]=v->getCoord(j);
		//printf("p[%ld]->(%f %f)\n",id,p.xyz[id+0],p.xyz[id+1]);
	      }
	  }  
	
	return true;
      } 
    
    p.xadj=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*(nRootNodes+1));
    p.vsize=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*(nRootNodes));
    p.vwgt=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*(nRootNodes));
    p.adjncy=(PartitionerIndex*)malloc(sizeof(PartitionerIndex)*nRootNodes*NNEI);

    p.xadj[0]=0;
    
    
    long long checkSum=0;
#pragma omp parallel for reduction(+:checkSum)
    for (long i=0;i<glb::num_omp_threads;i++)
      {	
	const network_iterator it_end=networkEnd();
	for (network_iterator it=networkBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {	  
	    const LocalIndex id = it->getLocalIndex()+1;
	    if (it->neighbors[NNEI-1]==NULL)
	      {
		p.xadj[id]=0;
		for (int j=0;j<NNEI-1;j++) 
		  if (it->neighbors[j]!=NULL) p.xadj[id]++;
	      }
	    else p.xadj[id]=NNEI;
	    // Redistribution cost (memory size)
	    p.vsize[id-1]= it->weight; 
	    // Computational cost
	    p.vwgt[id-1] = it->weight * pwFactor;
	    checkSum+=p.vwgt[id-1];
	  }
      }

    if ( mpiCom->sum(checkSum) > std::numeric_limits<PartitionerIndex>::max() )
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Graph weight sum exceeds the capacity of 'PartitionerIndex'.");
	glb::console->print<LOG_ERROR>("Try configuring parMetis to use 64bits index types.");
	glb::console->print<LOG_ERROR>("=> set IDXTYPEWIDTH to 64 in file metis.h and recompile ...");
	exit(-1);
      }

    for (unsigned long i=2;i<=nRootNodes;i++) p.xadj[i]+=p.xadj[i-1];

#pragma omp parallel for
    for (long i=0;i<glb::num_omp_threads;i++)
      {	
	const network_iterator it_end=networkEnd();
	for (network_iterator it=networkBegin(i,glb::num_omp_threads);it!=it_end;++it)
	  {	  
	    const LocalIndex id = it->getLocalIndex();
	    for (long j=p.xadj[id];j<p.xadj[id+1];j++)
	      {
		p.adjncy[j] = globalIdentityToGlobalIndex
		  (it->neighbors[j-p.xadj[id]]->getGlobalIdentity());
	      }
	  }
      }

    return true;
  }
  
  // repartition the tree
  // FIXME : post an Irecv before Isend and use waitall ...
  bool repart(std::vector<PartitionerIndex> &partition,       
	      RefinePartitionType type,
	      double tolerance,
	      double weightPerCell,
	      MpiCellDataExchangeT<Element,AnyNodeBase> &leavesExchange,
	      int nThreads=glb::num_omp_threads)	
  {
    const int myRank = mpiCom->rank();
    const int nParts = mpiCom->size();    

    ParmetisParams p(mpiCom);
    // drawGraph("GRAPH-PRE");
    glb::console->printFlush<LOG_INFO>("Repartitioning root nodes (%s) ... ",RefinePartitionTypeSelect().getString(type).c_str());

    glb::console->printFlush<LOG_PEDANTIC>("(graph) ");    
    if (!generateParmetisGraph(p,type,tolerance,weightPerCell)) return false;
    
    glb::console->printFlush<LOG_PEDANTIC>("(metis) ");    
    Partitioner::repart(p,type,partition);  

    /*
    PRINT_SRC_INFO(LOG_STD_ALL);    
    glb::console->print<LOG_STD_ALL>("THIS IS FOR DEBUG ONLY: DELETE ME !\n");
    {
      const unsigned long nRootNodes=getNRootNodes();
      partition.resize(nRootNodes);
      int nMin = glb::pParser->get<>("nReloc",
				     ParamsParser::defaultCategory(),
				     10,
				     "Number of cells to relocate per MPI process");
      for (network_iterator it=networkBegin();it!=networkEnd();++it)
	{
	  //partition[it->getLocalIndex()] = rand()%mpiCom->size();
	  
	  if (it->getLocalIndex()<nMin)
	    partition[it->getLocalIndex()] = (myRank+1)%mpiCom->size();//rand()%mpiCom->size();
	  else
	    partition[it->getLocalIndex()] = myRank;
	  
	}
    }
    */

    glb::console->printFlush<LOG_INFO>("done.\n");

    glb::console->printFlush<LOG_INFO>("Rebuilding tree structure ... ",RefinePartitionTypeSelect().getString(type).c_str());
    //glb::console->printFlush<LOG_PEDANTIC>("(tree) ");    

    std::vector<Root*> rootArr = getRootNodesArray();
    MpiCellDataExchangeT<Root,Root> rootExchange(mpiCom);
    std::vector< std::vector< MpiExchg_SharedTreeRoot > > toSend;
    unsigned long nSend=0;
    std::vector<unsigned long> nSendVec(nParts,0);
    std::vector< std::vector< typename MpiExchg_SharedTreeRoot::MpiStruct > > toReceive(nParts);
    unsigned long nReceive=0;    
    std::vector< std::vector<MPI_Request> > requests(4);    
        
    // DEBUG ONLY
    /*
    glb::console->print<LOG_STD_ALL>("REMOVE ME DEBUG IN SHAREDTREE.HXX");
    for (unsigned long i=0;i<partition.size();i++)
      partition[i]=myRank;

    if (myRank==0) {partition.back()=1;partition[11]=1;partition[12]=1;}
    if (myRank==1) {partition[10]=0;partition.back()=0;}
    */

    // compute how many roots will be sent to each process and in total   
    for (unsigned long i=0;i<partition.size();i++)
      if (partition[i]!=myRank)	  
	nSendVec[partition[i]]++;
      
    // setup rootExchange and reserve space to copy sent roots
    // ND: we must NOT underestimate the reserved size of sentRootArr[i]
    // as it would invalidate all pointers on resize !
    std::vector< std::vector<AnyNodeBase*> > sentRootTreeArr(nParts);    
    std::vector< std::vector<Root*> > sentRootPtrArr(nParts);    
    for (long i=0;i<nParts;i++)
      if (nSendVec[i]>0)
	{
	  rootExchange.sendRank.push_back(i);
	  sentRootTreeArr[i].reserve(nSendVec[i]);	  
	  sentRootPtrArr[i].reserve(nSendVec[i]);
	  nSend+=nSendVec[i];
	}
   
    // backup sent roots' trees for later use.
    // We cannot recycle them yet as this would invalidate neighbors
    // and we will still need them after the data is actually sent.
    std::vector<GlobalIdentity> freedGlobalIdentities;
    freedGlobalIdentities.reserve(nSend);
    for (unsigned long i=0;i<partition.size();i++)
      if (partition[i]!=myRank)
	{	
	  sentRootTreeArr[partition[i]].push_back(rootArr[i]->getChild());
	  sentRootPtrArr[partition[i]].push_back(rootArr[i]);
	  freedGlobalIdentities.push_back(rootArr[i]->getGlobalIdentity());
	  //sentRootHash.insert(std::make_pair(rootArr[i]->getGlobalIdentity().get(),&sentRootArr[partition[i]].back()));
	  //rootHash.erase(rootArr[i]->getGlobalIdentity().get());
	  rootArr[i]=NULL;
	  //rootPool.recycle(rootArr[i]);
	}  
    
    // we will pop_back the identities when reassigning them, and this should 
    // be done from lowest to highest !
    std::reverse(freedGlobalIdentities.begin(),freedGlobalIdentities.end());
    
    // untie the roots that will be sent from those that will stay
    // by replacing remaining roots displaced neighbors by shadows
    for (unsigned long i=0;i<rootExchange.sendRank.size();i++)
      {
	long index=rootExchange.sendRank[i];
	for (unsigned long j=0;j<sentRootPtrArr[index].size();j++)
	  {
	    Root* cur = sentRootPtrArr[index][j]; 
	    for (int k=0;k<NNEI;k++)
	      {
		Root* nei = cur->neighbors[k];
		if (nei==NULL) continue;
	      		
		//DenseHash_it it = rootHash.find(nei->getGlobalIdentity().get());
		//if (it!=rootHash.end())
		// check that the neighbor is local
		if (nei->getGlobalIdentity().rank() != myRank) continue;
		// and that it won't be sent to another process
		if (rootArr[nei->getGlobalIdentity().id()]!=NULL)
		  {
		    // this neighbor will stay on the current process
		    for (int l=0;l<NNEI;l++) 
		      {
			// find which of the neighbor's neighbors is the current 
			// sent root (cur)
			if (nei->neighbors[l]==cur)
			  {			    
			    // and replace it by a shadow ...
			    Root *sharedRoot;
			    shadowRootPool.pop(&sharedRoot);
			    sharedRoot->setGlobalIdentity(cur->getGlobalIdentity());
			    //sharedRoot->id = cur->id;
			    nei->neighbors[l]=sharedRoot;		
			    break;
			  }
		      }
		  }
	      }
	  }
      }

    // setup roots data to send
    toSend.resize(rootExchange.sendRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<toSend.size();i++)
      {	
	std::vector<Root*> &curSentRoot = sentRootPtrArr[rootExchange.sendRank[i]];
	toSend[i].resize(curSentRoot.size());
	for (unsigned long j=0;j<toSend[i].size();j++)
	  toSend[i][j].set(curSentRoot[j]);
      }

    // exchange information about which processes is sending to which process and how much
    const unsigned long nMax=mpiCom->max(rootExchange.sendRank.size());
    std::vector<long> sendCountInfo(nParts*nMax*2,-1);    
    for (unsigned long i=0;i<rootExchange.sendRank.size();i++)
      {
	sendCountInfo[myRank*nMax*2 + 2*i]=toSend[i].size();	
	sendCountInfo[myRank*nMax*2 + 2*i+1]=rootExchange.sendRank[i];
      }
    
    mpiCom->Allreduce_inplace(sendCountInfo,MPI_MAX);
    
    // Resize toReceive to hold the received roots
    unsigned long nTransferedRoots=0;
    unsigned long n=0;
    for (long i=0;i<nParts;++i)
      for (unsigned long j=0;j<nMax;j+=1,n+=2)
	{
	  if (sendCountInfo[n+1]>=0)	
	    nTransferedRoots += sendCountInfo[n];
	  if (sendCountInfo[n+1]==myRank)	    
	      toReceive[i].resize(sendCountInfo[n]);	  	    
	}
    
    // check from which ranks we are receiving and reorganize toReceive
    // so that it doe not have unnecessary entries
    for (long i=0;i<nParts;i++)
      if (toReceive[i].size()>0)
	{
	  nReceive+=toReceive[i].size();
	  std::swap(toReceive[i],toReceive[rootExchange.receiveRank.size()]);
	  rootExchange.receiveRank.push_back(i);	  
	}
    toReceive.resize(rootExchange.receiveRank.size());
  
    // Send info about the relocated sharedTree root nodes to the adequate processes
    rootExchange.exchangeStruct(toSend,toReceive,
				mpiType_sharedTreeRoot.getType(),false,
				mpiTagsStart_repart+0);    
    
    // Recycle the sent roots, we don't need them anymore
    // and we want to use their memory slots for the received ones
    // We keep their global identities in oldSentGid for later use

    std::vector< std::vector<GlobalIdentityValue> > oldSentGid(toSend.size());
    for (unsigned long i=0;i<toSend.size();i++)
      {
	long index=rootExchange.sendRank[i];
	oldSentGid[i].reserve(toSend.size());
	for (unsigned long j=0;j<toSend[i].size();j++)
	  {
	    //glb::console->print<LOG_PEDANTIC_ALL>("Sending %d.\n",sentRootPtrArr[index][j]->weight);
	    oldSentGid[i].push_back(sentRootPtrArr[index][j]->getGlobalIdentity().get());
	    rootPool.recycle(sentRootPtrArr[index][j]);
	  }
      } 
    sentRootPtrArr.clear();
   
    // Now create the received roots ,store them in a temporary array (receivedRootPtrArr)
    // and add their reference to the receivedRootHash table.
    std::vector< std::vector<Root*> > receivedRootPtrArr(toReceive.size());
    DenseHash receivedRootHash;    
    set_hash_empty_key(receivedRootHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(receivedRootHash,GlobalIdentity::max.get());
    const long receivedRootHashLoadGuess = nReceive;
    receivedRootHash.rehash(receivedRootHashLoadGuess/
			    (receivedRootHash.max_load_factor()-0.01));
    // glb::console->printToBuffer<LOG_PEDANTIC_ALL>("\n");
    for (unsigned long i=0;i<toReceive.size();i++)
      {
	receivedRootPtrArr[i].reserve(toReceive[i].size());	
	for (unsigned long j=0;j<toReceive[i].size();j++)
	  {
	    Root *root;
	    rootPool.pop(&root);
	 
	    if (freedGlobalIdentities.size()>0)
	      {
		root->setGlobalIdentity(freedGlobalIdentities.back());
		freedGlobalIdentities.pop_back();
	      }
	    else root->setGlobalIdentity(GlobalIdentity(myRank,rootPool.getUsedCount()-1));
	    // glb::console->printToBuffer<LOG_PEDANTIC_ALL>("+(%ld,%ld)",(long)root->id.rank(),(long)root->id.id());
	    root->weight = toReceive[i][j].weight;
	    //rootHash.insert(std::make_pair(root->id.get(),root));
	    receivedRootHash.insert(std::make_pair(toReceive[i][j].globalIdentity,root));
	    receivedRootPtrArr[i].push_back(root);
	  }
      }

    // glb::console->printToBuffer<LOG_PEDANTIC_ALL>("\n");
	    
    // localindex must be strictly growing consecutive integers starting at 0 at all time,
    // so if we received less roots than were sent, we have to change the global identity 
    // of the roots whose local index is higher than the number of local roots ...
    if (freedGlobalIdentities.size())
      {
	//local index must be lower than that !
	unsigned long nUsed=rootPool.getUsedCount(); 
	
	std::vector<GlobalIdentity> reassign;
	reassign.reserve(freedGlobalIdentities.size());
	for (unsigned long i=0;i<freedGlobalIdentities.size();i++)
	  if (freedGlobalIdentities[i].id()<nUsed) 
	    reassign.push_back(freedGlobalIdentities[i]);

	long curId=rootArr.size()-1;
	for (unsigned long i=0;i<reassign.size();++i)
	  {	   
	    while (rootArr[curId]==NULL) curId--;
	    // glb::console->printToBuffer<LOG_PEDANTIC_ALL>("-[(%ld,%ld)->(%ld,%ld)]",
	    // 					    (long)rootArr[curId]->id.rank(),(long)rootArr[curId]->id.id(),
	    // 					    (long)reassign[i].rank(),(long)reassign[i].id());
	    rootArr[curId--]->setGlobalIdentity(reassign[i]);
	    
	  }	
	// by construction, if (reassign.size()>0) then (curId < nUsed) here !
      }
    // glb::console->printToBuffer<LOG_PEDANTIC_ALL>("\n");
    // glb::console->flushBuffer<LOG_PEDANTIC_ALL>();

    // Now we can reconnect the received roots by correctly setting their neighbors
    for (unsigned long i=0;i<toReceive.size();i++)
      {
	for (unsigned long j=0;j<toReceive[i].size();j++)
	  {
	    Root *curRoot = receivedRootPtrArr[i][j];
	    // check each neighbors of curRoot
	    for (int k=0;k<NNEI;k++)
	      {
		if (toReceive[i][j].neighbors[k] == GlobalIdentity::empty.get()) 
		  {
		    curRoot->neighbors[k]=NULL;
		    continue;
		  }
		// REMEMBER: 'it' will be invalidated if we insert in rootHash !!!
		// this is due to dense hash maps implementation ...
		DenseHash_it it = receivedRootHash.find(toReceive[i][j].neighbors[k]);
		if (it == receivedRootHash.end())
		  {
		    // this neighbor was not received, it may be already local though.
		    //it = rootHash.find(toReceive[i][j].neighbors[k]);
		    //if (it == rootHash.end())		    
		    GlobalIdentity gid(toReceive[i][j].neighbors[k]);		    		   
		    if (gid.rank() != myRank)
		      {
			// No, it's not local (not received and from a different process)
			// -> this is a new shadow/ghost ...			
			// NB: they are duplicated and that's OK
			Root *sharedRoot;
			//#pragma omp critical 
			shadowRootPool.pop(&sharedRoot);
			// NOTE: we'll have to ask process gid.rank() for the
			// correct id later, this one is deprecated ...
			sharedRoot->setGlobalIdentity(toReceive[i][j].neighbors[k]); 
			curRoot->neighbors[k]=sharedRoot;
		      }
		    else if (rootArr[gid.id()]==NULL)
		      {
			// it was local but the neighbor root was sent to another node
			Root *sharedRoot;
			//#pragma omp critical 
			shadowRootPool.pop(&sharedRoot);
			// NOTE: we'll have to retrieve the
			// correct id later, this one is deprecated ...
			sharedRoot->setGlobalIdentity(toReceive[i][j].neighbors[k]); 
			curRoot->neighbors[k]=sharedRoot;
		      }
		    else
		      {
			// Yes, this neighbor is already local
			// lets tie curRoot to its new neighbor
			// -> rooArr index corresponds to the old local index
			//Root *nei = it->second;
			Root *nei = rootArr[gid.id()];
			curRoot->neighbors[k]=nei;
			// and also set curRoot as nei's neighbor
			for (int l=0;l<NNEI;l++)
			  {
			    // retieve which of the neighbors is curRoot
			    if (nei->neighbors[l]==NULL) continue;
			    // it used to be on a remote process ...
			    GlobalIdentity gid(nei->neighbors[l]->getGlobalIdentity());
			    if (gid.rank()==static_cast<unsigned int>(myRank)) 
			      continue;
			    it = receivedRootHash.find(gid.get());
			    if ((it!=receivedRootHash.end())&&
				(it->second == curRoot))
			      {
				shadowRootPool.recycle(nei->neighbors[l]);
				nei->neighbors[l]=curRoot;
				break;
			      }
			  }
		      }
		  }
		else 
		  {
		    // this is a newly received root, tie it to its neighbor
		    curRoot->neighbors[k] = it->second;
		    // don't need to set it->second's new neighbors, it will do it itself
		  }
	      }
	  }
      }  
    
    // We need to retrieve the new global ID's of the sent roots to update local 
    // shadows global ids and also be able to answer neighbors requests for new 
    // global ids (they may ask for a root that was sent).
    // NB: we sent to processes from which we received and receive from where we sent !
    //std::vector<MPI_Request> requests;    
    requests[0].resize(toReceive.size());
    std::vector< std::vector<GlobalIdentityValue> > newReceivedGid(toReceive.size());
    for (unsigned long i=0;i<toReceive.size();++i)
      {
	newReceivedGid[i].reserve(toReceive[i].size());
	for (unsigned long j=0;j<toReceive[i].size();++j)
	  newReceivedGid[i].push_back(receivedRootPtrArr[i][j]->getGlobalIdentity().get());
	  	
	mpiCom->Isend(&newReceivedGid[i][0],newReceivedGid[i].size(),
		      rootExchange.receiveRank[i],&requests[0][i],
		      mpiTagsStart_repart+1);	
      }

    // newSentGid is the new globalIdentity of sent roots
    // -> we can build a hash table mapping old GID to new GID ...
    std::vector< std::vector<GlobalIdentityValue> > newSentGid(toSend.size());
    DenseHashGidGid sentGidHash;
    set_hash_empty_key(sentGidHash,GlobalIdentity::empty.get());
    set_hash_deleted_key(sentGidHash,GlobalIdentity::max.get());
    const long sentGidHashLoadGuess = nSend;
    sentGidHash.rehash(sentGidHashLoadGuess/(sentGidHash.max_load_factor()-0.01));

    std::vector<int> indexOf(nParts,-1);
    for (unsigned long i=0;i<rootExchange.sendRank.size();++i)
      indexOf[rootExchange.sendRank[i]]=i;
    
    for (unsigned long i=0;i<toSend.size();++i)
      {
	// we already know how much we receive here ...
	int source=mpiCom->Probe(mpiTagsStart_repart+1); 
	int index = indexOf[source];
	newSentGid[index].resize(toSend[index].size());
	mpiCom->Recv(&newSentGid[index][0],toSend[index].size(),source,
		     mpiTagsStart_repart+1);
	for (unsigned long j=0;j<newSentGid[index].size();++j)
	  sentGidHash.insert(std::make_pair(oldSentGid[index][j],newSentGid[index][j]));
      }
    
    // Need to synchronize here as we don't want to risk messing up the following Isends 
    // with previous Probes ...  
    // --> actually, NOT NEEDED : we used MPI tags instead ...
    //mpiCom->barrier();
    //newReceivedGid.clear(); // NB: have to wait for the Isend to actually complete.
    
    // Finally, we also need to retrieve the remote neighbors' new global identity
    // -> first we have to identify the new neighborhood
    std::vector< std::vector<Root*> > shadowToReceive(nParts);
    std::vector< long > shadowReceiveRank;    
    network_iterator it_end=shadowNetworkEnd();
    for (network_iterator it=shadowNetworkBegin();it!=it_end;++it)	
      {
	GlobalIdentity gid(it->getGlobalIdentity());
	if (gid.rank() != myRank) 
	  shadowToReceive[gid.rank()].push_back(*it);	 
      }

    for (long i=0;i<nParts;i++)
      {
	if (shadowToReceive[i].size()>0) 
	  shadowReceiveRank.push_back(i); 	   
      }
    
    const long nShadowMax=mpiCom->max(shadowReceiveRank.size());
    std::vector<long> shadowReceiveCountInfo(nParts*nShadowMax*2,-1);    
    for (unsigned long i=0;i<shadowReceiveRank.size();i++)
      {
	long index = shadowReceiveRank[i];
	shadowReceiveCountInfo[myRank*nShadowMax*2 + 2*i]=shadowToReceive[index].size();
	shadowReceiveCountInfo[myRank*nShadowMax*2 + 2*i+1]=index;
      }

    // now everybody should know how much to send, as well as the destination
    mpiCom->Allreduce_inplace(shadowReceiveCountInfo,MPI_MAX);
    
    std::vector< std::vector<GlobalIdentityValue> > shadowGidToSend(nParts);
    std::vector< long > shadowSendRank;
    n=0;    
    for (long i=0;i<nParts;++i)
      for (long j=0;j<nShadowMax;j+=1,n+=2)
	if (shadowReceiveCountInfo[n+1]==myRank)
	  {
	    shadowGidToSend[i].resize(shadowReceiveCountInfo[n]);	
	    shadowSendRank.push_back(i);
	  }

    // Now that they are allocated, prepare to receive the shadowGid to send
    std::vector<MPI_Request> shadowGidToSendRequests(shadowSendRank.size());
    for (unsigned long i=0;i<shadowSendRank.size();i++)
      {
	int source = shadowSendRank[i];
	mpiCom->Irecv(&shadowGidToSend[source][0],shadowGidToSend[source].size(),source,
		      &shadowGidToSendRequests[i],mpiTagsStart_repart+2);
      }
    
    // first send the list of roots we want to receive the globalidentity of
    std::vector< std::vector<GlobalIdentityValue> > shadowGidToReceive(nParts);
    requests[1].resize(shadowReceiveRank.size());
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<shadowReceiveRank.size();i++)
      {
	long index = shadowReceiveRank[i];
	shadowGidToReceive[index].reserve(shadowToReceive[index].size());
	for (unsigned long j=0;j<shadowToReceive[index].size();j++)
	  shadowGidToReceive[index].push_back(shadowToReceive[index][j]->getGlobalIdentity().get());

	mpiCom->Isend(&shadowGidToReceive[index][0],shadowGidToReceive[index].size(),
		      index,&requests[1][i],mpiTagsStart_repart+2);
      }

    // Use latency to update the shadow root neighbors that were local but are not anymore 
    // --> sent roots
    for (network_iterator it=shadowNetworkBegin();it!=it_end;++it)	
      {
	GlobalIdentity gid(it->getGlobalIdentity());
	if (gid.rank() == myRank) 
	  {
	    DenseHashGidGid_it it2 = sentGidHash.find(gid.get());
	    if (it2!=sentGidHash.end())
	      it->setGlobalIdentity(it2->second);
	  }
      }
    
    // wait to receive the shadowGidToSend
    mpiCom->Waitall(shadowGidToSendRequests);
    // wait for data in shadowGidToReceive to be sent, so that we can reuse the buffer !
    mpiCom->Waitall(requests[1]); 

    // prepare to receive
    std::vector<MPI_Request> shadowGidToReceiveRequests(shadowReceiveRank.size());
    for (unsigned long i=0;i<shadowReceiveRank.size();i++)
      {
	int source=shadowReceiveRank[i];
	mpiCom->Irecv(&shadowGidToReceive[source][0],shadowGidToReceive[source].size(),
		      source,&shadowGidToReceiveRequests[i],mpiTagsStart_repart+3);
      }
     
    // then look-up the requested roots and send back their new global identity
    requests[2].resize(shadowSendRank.size());
    for (unsigned long i=0;i<shadowSendRank.size();i++)
      {
	//int source=mpiCom->Probe(1);
	//mpiCom->Recv(&shadowGidToSend[source][0],shadowGidToSend[source].size(),source,1);
	int source = shadowSendRank[i];
	for (unsigned long j=0;j<shadowGidToSend[source].size();++j)
	  {
	    long id=GlobalIdentity(shadowGidToSend[source][j]).id();
	    if (rootArr[id]!=NULL)
	      shadowGidToSend[source][j] = rootArr[id]->getGlobalIdentity().get();
	    else
	      shadowGidToSend[source][j] = sentGidHash.find(shadowGidToSend[source][j])->second;	      
	  }
	mpiCom->Isend(&shadowGidToSend[source][0],shadowGidToSend[source].size(),source,
		      &requests[2][i],mpiTagsStart_repart+3);
      }

    // wait to receive the shadowGidToSend
    mpiCom->Waitall(shadowGidToReceiveRequests);
    // receive the new global identities and assign them to the corresponding shadow roots
    // after that we are done with updating the global identities
    for (unsigned long i=0;i<shadowReceiveRank.size();i++)
      {
	int source = shadowReceiveRank[i];
	//int source=mpiCom->Probe(2);	
	//mpiCom->Recv(&shadowGidToReceive[source][0],shadowGidToReceive[source].size(),source,2);
	for (unsigned long j=0;j<shadowGidToReceive[source].size();++j)
	  shadowToReceive[source][j]->setGlobalIdentity(shadowGidToReceive[source][j]);
      }    
    
    // Finally, retrieve and exchange the tree structures 
    // compactTree *MUST* be of an integral type big enough to store node data in it
    // long is actually bigger than necessary, so use int (short would probably work ...)
    typedef unsigned int CompactTreeData;
    std::vector< std::vector<CompactTreeData> > compactTree(toSend.size());     
    requests[3].resize(toSend.size());
    // compute how many nodes/leaves in the tree, build a compact tree data and Isend it
#pragma omp parallel for num_threads(nThreads)
    for (unsigned long i=0;i<toSend.size();i++)
      {
	long index=rootExchange.sendRank[i];
	char treeIndex=0;
	long size=0;
	std::vector<CompactTreeData> nodeData;
	nodeData.reserve(sentRootTreeArr[index].size()*2);
	for (unsigned long j=0;j<sentRootTreeArr[index].size();j++)
	  size+=sentRootTreeArr[index][j]->
	    buildCompactTree(compactTree[i],treeIndex,nodeData);

	// Add node data at the end of the compact tree so that we can transfer everything
	// in a single pass.
	unsigned long cTSize=compactTree[i].size();
	compactTree[i].reserve(compactTree[i].size()+nodeData.size()+1);
	compactTree[i].insert(compactTree[i].end(),nodeData.begin(),nodeData.end());
	compactTree[i].push_back(cTSize);

	mpiCom->Isend(&compactTree[i][0],compactTree[i].size(),
		      rootExchange.sendRank[i],&requests[3][i],mpiTagsStart_repart+4);
      }  

    //std::vector< std::vector<Element*> > sentLeafArr(nParts);
    //sentLeafArr.resize(nParts);
    leavesExchange.send.clear();
    leavesExchange.send.resize(nParts);
    leavesExchange.sendRank.clear();
    // lets free the sent nodes and store the leaves to send while 
    // compact trees are being sent so that latency is covered and
    // we can reuse recycled nodes ...
    for (unsigned long i=0;i<rootExchange.sendRank.size();i++)
      {
	const long index=rootExchange.sendRank[i];
	for (unsigned long j=0;j<sentRootTreeArr[index].size();j++)
	  {
	    Node *node=sentRootTreeArr[index][j]->
	      getLeavesRecycleNodes(leavesExchange.send[index],nodePool);
		    
	    if (node!=NULL) nodePool.recycle(node);
	  }
	leavesExchange.sendRank.push_back(index);
      }
    sentRootTreeArr.clear();

    // and now receive the trees and rebuild them on the local node
    std::vector< std::vector<CompactTreeData> > receivedCompactTree(toReceive.size());
    //std::vector< std::vector<Leaf*> > receivedLeafParentsArr(toReceive.size());
    //receivedLeafParentsArr.resize(toReceive.size());
    leavesExchange.receive.resize(nParts);
    leavesExchange.receiveRank.clear();

    indexOf.assign(indexOf.size(),-1);
    for (unsigned long i=0;i<rootExchange.receiveRank.size();++i)
      indexOf[rootExchange.receiveRank[i]]=i;
    for (unsigned long i=0;i<toReceive.size();i++)
      {
	int count=0;
	int source=mpiCom->
	  template ProbeCount<CompactTreeData>(count,mpiTagsStart_repart+4);
	int index = indexOf[source];

	receivedCompactTree[index].resize(count);
	mpiCom->Recv(&receivedCompactTree[index][0],count,source,mpiTagsStart_repart+4);
	restoreCompactTree(receivedRootPtrArr[index],
			   receivedCompactTree[index],
			   receivedCompactTree[index].begin()+
			   receivedCompactTree[index].back(),
			   leavesExchange.receive[source]);
	
	leavesExchange.receiveRank.push_back(source);
      }

    // wait for all pending Isend requests: locally allocated buffers
    // may be freed before they complete when function returns ...
    for (unsigned long i=0;i<requests.size();i++)
      for (unsigned long j=0;j<requests[i].size();j++)
	mpiCom->Wait(&requests[i][j]);

    // drawGraph("GRAPH-POST");

    // glb::console->printFlush<LOG_PEDANTIC>("(graph2) ");    
    // if (!generateParmetisGraph(p,tolerance)) return false;
   
    // glb::console->printFlush<LOG_INFO>("(metis2) ");
    // partition = Partitioner::repart(p,type);  

    glb::console->print<LOG_INFO>("done. (%ld roots relocated)\n",nTransferedRoots);
    
    return true;
  }

  template <class IT>
  Node *restoreTree_recurse(const std::vector<IT> &compactTree,
			    typename std::vector<IT>::const_iterator &ctIt,
			    long &index,int &k,std::vector<AnyNodeBase *> &parents)
  {
    if (compactTree[index]&(1ul<<k))
      {
	if ((++k) == sizeof(IT)*8) {k=0;++index;}
	Node *node;
	nodePool.pop(&node);
	node->setChildrenIndex(*ctIt);++ctIt;
	node->children[0]=restoreTree_recurse(compactTree,ctIt,index,k,parents);	
	if (node->children[0] == NULL) parents.push_back(node);
	else node->children[0]->parent = node;
	node->children[1]=restoreTree_recurse(compactTree,ctIt,index,k,parents);
	if (node->children[1] == NULL) parents.push_back(node);
	else node->children[1]->parent = node;
	return node;
      }
    else
      {	
	if ((++k) == sizeof(IT)*8) {k=0;++index;}
	return NULL;
      }
    return NULL;
  }

  template <class IT>
  void restoreCompactTree(std::vector<Root*> &newRoots, 
			  const std::vector<IT> &compactTree, 
			  typename std::vector<IT>::const_iterator ctIt,
			  std::vector<AnyNodeBase *> &parents)
  {
    long index=0;
    int k=0;
    typename std::vector<IT>::const_iterator it = ctIt;
    for (unsigned long i=0;i<newRoots.size();++i)
      {
	newRoots[i]->setChild(restoreTree_recurse(compactTree,it,index,k,parents));
	if (newRoots[i]->getChild()==NULL) parents.push_back(newRoots[i]);
	else newRoots[i]->getChild()->setParent(newRoots[i]);
      }
  }

  long countCompactTreeLeaves(const std::vector<unsigned long> &compactTree)
  {
    long size=0;
    for (unsigned long i=0;i<compactTree.size();i++)
      for (size_t j=0;j<sizeof(unsigned long)*8;j++)
	if (compactTree[i]&(1ul<<j)) ++size;
    return size;
  }

  void drawGraph(const char *fname)
  {
    const int myRank = mpiCom->rank();
    const int nParts = mpiCom->size();
    updateRootNodesCount();
    glb::console->print<LOG_STD>("\n");
    glb::console->printFlush<LOG_STD>("Drawing the graph ... ");
    char newname[255];
    sprintf(newname,"%s.diag.dot",fname);
    FILE *f;
    if (myRank == 0)
      {
	f=fopen(newname,"w");
	fprintf(f,"digraph G {");
      }
    else
      {
	int i=0;
	mpiCom->Recv(&i,1,myRank-1,mpiTagsStart_general+0);
	f=fopen(newname,"r+");
	fseek(f,0,SEEK_END);
      }
    bool rankPrint=glb::console->getRankPrint();
    glb::console->setRankPrint(false);
    glb::console->setStartingLine(false);
    glb::console->printFlush<LOG_STD_ALL>("(%d) ",myRank);    

    //int i=0;
    //int j=0;
    
    std::vector<std::string> boxtype;
    boxtype.push_back(std::string("box"));
    boxtype.push_back(std::string("circle"));
    boxtype.push_back(std::string("diamond"));
    boxtype.push_back(std::string("triangle"));

    std::vector<std::string> color;
    color.push_back(std::string("blue"));
    color.push_back(std::string("green"));
    color.push_back(std::string("orange"));
    color.push_back(std::string("red"));

    std::vector<std::string> style;
    style.push_back(std::string("dashed"));
    style.push_back(std::string("diagonals"));

    // printf("START\n");
    // printf("NUSED : %ld\n",rootPool.getUsedCount());
    const network_iterator it_end=networkEnd();
    for (network_iterator it=networkBegin();it!=it_end;++it)	   
      {
	int boundary=0;
	int shared=0;
	for (int i=0;i<NNEI;i++)
	  {
	    if (it->neighbors[i]==NULL) boundary=1;
	    else if (it->neighbors[i]->getGlobalIdentity().rank()!=myRank) shared=1;
	  }
	
	// printf("GID : %ld %ld\n",(long)it->getGlobalIdentity().rank(),(long)it->getGlobalIdentity().id());
	GlobalIndex a=globalIdentityToGlobalIndex(it->getGlobalIdentity());
	// printf("->%ld\n",a);
	/*
	printf("%ld [shape=%s,color=%s,style=%s];\n",a,
	       boxtype[boundary].c_str(),
	       color[myRank].c_str(),
	       style[shared].c_str());
	*/
	
	fprintf(f,"%ld [shape=%s,color=%s,style=%s];\n",a,
		boxtype[boundary].c_str(),
		color[myRank].c_str(),
		style[shared].c_str());

	
	//nodeID.insert(std::pair<double,int>(n->uniqueID(),i));
      }

    //printf("STOP\n");
    //std::vector<node_it> nei;
    //std::vector<node_it> id;
    //bool linked;
    char txt[1024];
    //std::set<node_it,NDcomplex_node::compareItLess> uniq;
    //std::vector<node_it>::iterator nei_it;
    // i=0;
    
    for (network_iterator it=networkBegin();it!=it_end;++it)	   
      {
	for (int i=0;i<NNEI;i++)
	  {
	    if (it->neighbors[i] == NULL) continue;
	    strcpy(txt,"");	    
	    if (it->neighbors[i]->getGlobalIdentity().rank()!=myRank)
	      sprintf(txt,"style=dotted,color=%s",color[myRank].c_str());
	    else
	      sprintf(txt,"color=%s",color[myRank].c_str());

	    GlobalIndex a=globalIdentityToGlobalIndex(it->getGlobalIdentity());
	    GlobalIndex b=globalIdentityToGlobalIndex(it->neighbors[i]->getGlobalIdentity());

	    fprintf(f,"%ld->%ld [%s];\n",a,b,txt);
	  }
      }
    
    if (myRank == nParts-1)
      {
	fprintf(f,"}");
	fclose(f);      
	char cmd[255];
	//sprintf(cmd,"fdp -Tps %s.diag.dot -o %s.diag.ps; gv %s.diag.ps &",fname,fname,fname);
	sprintf(cmd,"fdp -Tps %s.diag.dot -o %s.diag.ps",fname,fname);
	int err = system(cmd);	
	if (err==-1) 
	  glb::console->print<LOG_WARNING>("Call to system command '%s' failed.\n",cmd);
	glb::console->print<LOG_STD_ALL>("done.\n");	
      }
    else
      {
	int i=1;
	fclose(f);      
	mpiCom->Send(&i,1,myRank+1,mpiTagsStart_general+0);
      }
    mpiCom->barrier();
    glb::console->setRankPrint(rankPrint);
    glb::console->setStartingLine(true);
  }

};

/** \}*/
#include "../internal/namespace.footer"
#endif
