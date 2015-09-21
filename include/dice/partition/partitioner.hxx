#ifndef __PARTITIONER_HXX__
#define __PARTITIONER_HXX__

#include <limits>
//#include <parmetis.h>

#include "../dice_globals.hxx"
#include "../partition/partitionType.hxx"
#include "../partition/parmetisParams.hxx"
#include "../geometry/boundaryType.hxx"


#include "../internal/namespace.header"

class Partitioner
{
public:
  typedef ParmetisParams::Index Index;
  typedef ParmetisParams::Float Float;

  template <class SGT>
  static long buildFromGrid(const SGT *sg, 
			    std::vector<Index> &partition, 
			    PartitionType type, 
			    double tolerance=1.05, 
			    MpiCommunication *com=glb::mpiComWorld)
  {
    typedef typename SGT::Cell SGCell;
    ParmetisParams p(com);
          
    if (com->size()==1)
      {
	partition.clear();
	return 0;
      }
    p.setType(type);    
    p.setTolerance(tolerance);
    p.nDims=sg->getNDims();   

    glb::console->printFlush<LOG_PEDANTIC>("(graph) ");
    std::vector<unsigned long> nCells=sg->getNCells();
    std::vector<Index> xadj;
    std::vector<Index> adjncy;
    std::vector<Index> vtxdist;
    std::vector<Float> xyz;
    
    if (nCells[p.nDims] > static_cast<unsigned long>(std::numeric_limits<Index>::max()))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("Integer type in ParMetis is too short (%ld bytes)\n",sizeof(Index));
	glb::console->print<LOG_ERROR>
	  ("See ParMetis documentation and change it to 8 bytes.\n",sizeof(Index));
	exit(-1);
      }
      
    vtxdist.resize(com->size()+1);    
    for (long i=0;i<=com->size();i++)
      {
	vtxdist[i]=(i*nCells[p.nDims])/com->size();	
      }
    vtxdist.back()=nCells[p.nDims];
    p.vtxdist = &vtxdist[0];

    const unsigned long i0=vtxdist[com->rank()];
    const unsigned long imax=vtxdist[com->rank()+1];  
    const unsigned long nSimplices=imax-i0;
    
    if (type != PartitionTypeV::PH)
      {
	xadj.resize(nSimplices + 1);
	adjncy.resize(xadj.size()*(p.nDims+1));
	
	xadj[0]=0;
	if (glb::num_omp_threads<3)
	  {
	    for (unsigned long i=i0;i<imax;i++)
	      {
		std::vector<SGCell> cellV;
		sg->getNeighbors(SGCell(p.nDims,i),cellV);//std::back_inserter(cellV));	
		xadj[i-i0+1]=xadj[i-i0]+cellV.size();	
		for (unsigned long j=0;j<cellV.size();j++) 
		  adjncy[xadj[i-i0]+j] = cellV[j].id();				 
	      }	   
	  }
	else
	  {
#pragma omp parallel for
	    for (unsigned long i=i0;i<imax;i++)
	      {
		std::vector<SGCell> cellV;
		sg->getNeighbors(SGCell(p.nDims,i),cellV);//std::back_inserter(cellV));	
		xadj[i-i0+1]=cellV.size();	
	      }

	    for (unsigned long i=2;i<=nSimplices;i++) xadj[i]+=xadj[i-1];

#pragma omp parallel for
	    for (unsigned long i=i0;i<imax;i++)
	      {
		std::vector<SGCell> cellV;
		sg->getNeighbors(SGCell(p.nDims,i),cellV);//std::back_inserter(cellV));	
		for (unsigned long j=0;j<cellV.size();j++) 
		  adjncy[xadj[i-i0]+j] = cellV[j].id();
	      }
	  }

	adjncy.resize(xadj.back()); 
	p.xadj=&xadj[0];
	p.adjncy=&adjncy[0];
      }   

    if (type != PartitionTypeV::KWAY)
      {
	xyz.resize(nCells[p.nDims]*p.nDims);
	for (unsigned long i=i0;i<imax;i++)
	  {
	    double tmp[p.nDims];
	    sg->getPosition(SGCell(p.nDims,i),&tmp[0]);
	    std::copy(tmp,tmp+p.nDims,&xyz[(i-i0)*p.nDims]);
	    //sg->getPosition(SGCell(p.nDims,i),&xyz[(i-i0)*p.nDims]);
	  }
	p.xyz=&xyz[0];
      }      

    glb::console->printFlush<LOG_PEDANTIC>
      ("(%s) ",PartitionTypeSelect().getString(p.type).c_str()); 
  
    //partition=build(p,false);
    p.partition(partition);
  
    xadj.clear();
    adjncy.clear();
    vtxdist.clear();
    xyz.clear();
    
    glb::console->printFlush<LOG_PEDANTIC>("(gathering) ");
    communicatePartition(partition,i0,com);
    
    return p.edgeCut;
  }  
  /*
  static std::vector<Index> build(ParmetisParams &p, bool communicate=true)
  { 
    
    std::vector<Index> partition;
    unsigned long nNodes=p.vtxdist[p.mpiCom->rank()+1]-p.vtxdist[p.mpiCom->rank()];
    MPI_Comm cm=p.mpiCom->getCom(); // communicator

    partition.resize(nNodes);

    if (p.type == PartitionTypeV::KWAY)
      {
	// Purely topological version
	// Slow but gives a good partition
	// Prefered version if you have enough memory

	ParMETIS_V3_PartKway(p.vtxdist,p.xadj,p.adjncy,
			     p.vwgt,p.adjwgt,&p.wgtflag,
			     &p.numflag,&p.ncon,&p.nParts,
			     p.tpwgts,p.ubvec,
			     p.opt,&p.edgeCut,&partition[0],&cm);

      }
    else if (p.type == PartitionTypeV::PH)
      {
	// Purely geometrical, using PH-curve
	// very fast and cheap 
	// Sometimes the partition is good, sometimes not ...

	ParMETIS_V3_PartGeom(p.vtxdist,&p.nDims,p.xyz,&partition[0],&cm);
	p.edgeCut=-1;
      }
	
    else if (p.type == PartitionTypeV::PH_KWAY)
      {
	// Hybrid method
	// Not sure if the quality is better or worse than KWAY
	// uses more memory, but faster
	   	
	ParMETIS_V3_PartGeomKway(p.vtxdist,p.xadj,p.adjncy,
				 p.vwgt,p.adjwgt,&p.wgtflag,
				 &p.numflag,&p.nDims,p.xyz,&p.ncon,&p.nParts,
				 p.tpwgts,p.ubvec,
				 p.opt,&p.edgeCut,&partition[0],&cm);
      }
    else
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Unknown partition type.\n");
	exit(-1);
      }
  
    if (communicate)
      communicatePartition(partition,p.vtxdist[p.mpiCom->rank()],p.mpiCom);   

    return partition;
  }
  */
  //static std::vector<Index> 
  static void repart(ParmetisParams &p, RefinePartitionType type, std::vector<Index> &partition, 
		     bool communicate=false)
  { 
    partition.clear();
    //std::vector<Index> partition;
    
    //unsigned long nNodes=p.vtxdist[p.mpiCom->rank()+1]-p.vtxdist[p.mpiCom->rank()];
    //MPI_Comm cm=p.mpiCom->getCom(); // communicator
    
    //partition.resize(nNodes);
    
    p.setRefineType(type);    
    /*
    std::vector<Float> xyz;
    if (type == RefinePartitionTypeV::PH)
      {
	xyz.resize(nCells[p.nDims]*p.nDims);
	for (unsigned long i=i0;i<imax;i++)
	  sg->getPosition(SGCell(p.nDims,i),&xyz[(i-i0)*p.nDims]);
	  
	p.xyz=&xyz[0];
      }      
    */

    p.refine(partition);

    //p.xyz=NULL;
    /*
    ParMETIS_V3_AdaptiveRepart(p.vtxdist,p.xadj,p.adjncy,
			       p.vwgt,p.adjwgt,p.vsize, &p.wgtflag,
			       &p.numflag,&p.ncon,&p.nParts,
			       p.tpwgts,p.ubvec,&p.itr,
			       p.opt,&p.edgeCut,&partition[0],&cm);
    */
    if (communicate)
      communicatePartition(partition,p.vtxdist[p.mpiCom->rank()],p.mpiCom);			   

    //return partition;
  }

  /*
  template <class LMT>
  static long improve(const LMT *localMesh, std::vector<Index> &partition, double tolerance=1.05, MpiCommunication *com=glb::mpiComWorld)
  {
    const Index nParts=com->size(); // number fo partitions    
    const Index ndims=LMT::NDIM;
    std::vector<unsigned long> nCells=localMesh->getNCells();
    const unsigned long nSimplices=nCells[ndims];
    std::vector<Index> xadj;
    std::vector<Index> adjncy;
    std::vector<Index> vtxdist;
    std::vector<Float> xyz;
    
    std::vector<unsigned long> i0V(nParts,0);
    std::vector<unsigned long> imaxV(nParts,0);
    imaxV[com->rank()]=nCells[ndims];
    com->Allreduce_inplace(imaxV,MPI_SUM); // to avoid buggy MPI_Allgather
    for (int i=1;i<nParts;i++) 
      {
	i0V[i]=imaxV[i-1];
	imaxV[i]+=imaxV[i-1];
      }
    const unsigned long i0=i0V[com->rank()];
    const unsigned long nSimplicesTotal=imaxV.back();
    //const unsigned long imax=imaxV[com->rank()]; 
   
    if (nSimplicesTotal > std::numeric_limits<Index>::max())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Integer type in ParMetis is too short (%ld bytes)\n",sizeof(Index));
	glb::console->print<LOG_ERROR>("See ParMetis documentation and change it to 8 bytes.\n",sizeof(Index));
	exit(-1);
      }
    

    vtxdist.resize(nParts+1,0);    
    for (unsigned long i=1;i<=nParts;i++)
      vtxdist[i]=imaxV[i-1];	
        
    Index wgtflag=0; // No weight
    Index numflag=0; // C-style
    Index ncon=1; // number of constraints    
    Float tpwgts[ncon*nParts]; // weight expected on each partition (1/nParts)
    Float ubvec[ncon]; // tolerance on repartition (1.05 is fine)
    Index opt[3]={1,0,0}; // options [use options,output timing,seed]        
    Index edgeCut; // used as output (number of cut edges)
    MPI_Comm cm=com->getCom(); // communicator

    glb::console->printFlush<LOG_INFO>("(graph) ");

    for (unsigned long j=0;j<nParts*ncon;j++)
      tpwgts[j]=1.0/nParts;

    for (unsigned long j=0;j<ncon;j++) ubvec[j]=tolerance;
      
    partition.resize(nSimplices);   

    xadj.resize(nSimplices + 1);
    adjncy.resize(xadj.size()*(ndims+1));
	
    xadj[0]=0;

    typedef typename LMT::simplex_iterator simplex_iterator;   
#pragma omp parallel for
    for (unsigned long th=0;th<glb::num_omp_threads;th++)
      {
	const simplex_iterator it_end=localMesh->simplexEnd();
	for (simplex_iterator it=localMesh->simplexBegin(th,glb::num_omp_threads);it!=it_end;it++)
	  {	
	    unsigned long curId=it->getLocalIndex();
	    int nnei=it->getNeighborsCount();
	    xadj[curId+1]=nnei;
	  }
      }

    for (unsigned long i=2;i<=nSimplices;i++) xadj[i]+=xadj[i-1];

#pragma omp parallel for
    for (unsigned long th=0;th<glb::num_omp_threads;th++)
      {
	const simplex_iterator it_end=localMesh->simplexEnd();
	for (simplex_iterator it=localMesh->simplexBegin(th,glb::num_omp_threads);it!=it_end;it++)
	  {
	    unsigned long neiId[LMT::Simplex::NNEI];
	    unsigned long curId=it->getLocalIndex();
	    int nnei=it->getNeighborsLocalIndex(neiId);
	    for (unsigned long j=0;j<nnei;j++) 
	      adjncy[xadj[curId]+j] = neiId[j];
	  }	
      }
    
    adjncy.resize(xadj.back()); 
             
    glb::console->printFlush<LOG_INFO>("(refine-KWAY) ");
    
    ParMETIS_V3_RefineKway(&vtxdist[0],&xadj[0],&adjncy[0],
			   NULL,NULL,&wgtflag,
			   &numflag,&ncon,&nParts,
			   tpwgts,ubvec,
			   opt,&edgeCut,&partition[0],&cm);    
    xadj.clear();
    adjncy.clear();
    vtxdist.clear();
    
    glb::console->printFlush<LOG_INFO>("(gathering) ");
    communicatePartition(partition,i0,com);
  }
  */

  // FIXME: this could be significantly improved when the number of processes is large ...
  // after calling this, partition will contain the indices belonging to the current node only
  static void communicatePartition(std::vector<Index> &partition, unsigned long i0, 
				   MpiCommunication *com=glb::mpiComWorld,
				   std::vector<Index> *ref=NULL)
  {
    unsigned long nPart = com->size();
    unsigned long rank = com->rank();

    //std::vector<Index> sendBuf[nPart];
    std::vector< std::vector<Index> > sendBuf(nPart);

    std::vector<unsigned long> count(nPart*nPart,0); // count[i*nPart+j] -> how many i should send to j
 
    if (ref==NULL)
      {
	for (unsigned long i=0;i<partition.size();i++)
	  {
	    count[rank*nPart+partition[i]]++;
	    sendBuf[partition[i]].push_back((Index)(i0+i));
	  } 
      }
    else
      {
	for (unsigned long i=0;i<partition.size();i++)
	  {
	    count[rank*nPart+partition[i]]++;
	    sendBuf[partition[i]].push_back((*ref)[i]);
	  } 
      }

  
   // NOTE: there is a bug in MPI_Allgather of openMPI v1.4.x
   // so we use allreduce here as performence is not crucial
   //com->Allgather_inplace(count);
   com->Allreduce_inplace(count,MPI_SUM);

   // this tests for the bug 
   /*
     glb::console->print<LOG_DEBUG>(" MPI_Allgather is %sbugged.\n",com->AllGatherIsBugged()?"":"NOT ");
     
     for (unsigned long i=0;i<nPart*nPart;i++) 
     {
     if (rank==(i/nPart)) count[i]=com->rank();
     else count[i]=0;
     }
   
     com->Allgather_inplace(count);
   
     for (int k=0;k<com->size();k++)
     {
     for (unsigned long i=0;i<nPart;i++) 
     for (unsigned long j=0;j<nPart;j++)
     glb::console->printToBuffer<LOG_DEBUG>("%ld,",count[i*nPart+j]);
     glb::console->printToBuffer<LOG_DEBUG>("]\n");
     glb::console->flushBuffer<LOG_DEBUG>();
     com->barrier();
     }
     //exit(0);
   */
   //glb::console->print<LOG_STD>("rank %d: reached barrier 1\n",com->rank());
   //com->barrier();
    
    unsigned long mySize=0;
    for (unsigned long i=0;i<nPart;i++)
      mySize+=count[i*nPart+rank];
    
    partition.resize(mySize); // receiveBuf 
    Index *rcv=&partition[0]; 
   
    for (unsigned long i=0;i<nPart;i++) // sender
      {		
	for (unsigned long j=0;j<nPart;j++) // receiver
	  {	
	    long nExchg=count[i*nPart+j];
	    if (nExchg==0) continue; // nobody does anything
	    if ((i!=rank)&&(j!=rank)) continue; // I do nothing ...
	    
	    if (i==j)
	      {
		Index *snd=&sendBuf[j][0];

		glb::console->printNewLine<LOG_DEBUG>
		  ("Sending %ld elements to myself\n",nExchg);
		std::copy(snd,snd+nExchg,rcv);
		//glb::console->printNewLine<LOG_DEBUG>("Sending %ld elements to myself DONE\n",nExchg);
		rcv+=nExchg;
	      }
	    else if (i==rank) // I send
	      {
		Index *snd=&sendBuf[j][0];

		glb::console->printNewLine<LOG_DEBUG>
		  ("Sending %ld element to %ld\n",nExchg,j);
		com->Send(snd,nExchg,j);
		//glb::console->printNewLine<LOG_DEBUG>("Sending %ld element to %ld DONE\n",nExchg,j);
	      }
  	    else // I receive	
	      {
		glb::console->printNewLine<LOG_DEBUG>("Receiving %ld element from %ld\n",nExchg,i);	
		com->Recv(rcv,nExchg,i);
		//glb::console->printNewLine<LOG_DEBUG>("Receiving %ld element from %ld DONE\n",nExchg,i);	
		rcv+=nExchg;	
	      }
	    
	  }		
      }  

   
    //glb::console->print<LOG_DEBUG>("rank %d: reached barrier\n",com->rank());
    //com->barrier();
 
  }

};

#include "../internal/namespace.footer"
#endif
