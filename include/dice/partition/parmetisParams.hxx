#ifndef __PARMETIS_PARAMS_HXX__
#define __PARMETIS_PARAMS_HXX__

#include <stdlib.h>
//#include <parmetis.h>

#include "../dice_globals.hxx"
#include "../partition/parmetisInterface.hxx"
#include "../partition/partitionType.hxx"


#include "../internal/namespace.header"

struct ParmetisParams
{
  typedef idx_t Index;
  typedef real_t Float;
  
  Index nParts; // Number of partitions
  Index ncon; // number of constraints 
  
  Index wgtflag; // Amout of work associated to a node  
  Index numflag; // C-style   
  Float *tpwgts; // weight expected on each partition (1/nParts)
  Float *ubvec; // tolerance on repartition (1.05 is fine)
  Index opt[4]; // options [use options,output timing,seed]   
  Index *vsize; // for adaptiverepat, amount of memory associated to a node
  Float itr;
  PartitionType type;
  RefinePartitionType refineType;
  MpiCommunication *mpiCom;

  Index *vwgt;
  Index *adjwgt;
  Index *xadj;
  Index *adjncy;
  Index *vtxdist; 
  Float *xyz; 
  Index nDims;

  bool owner;

  Index edgeCut; // output parameter
  
  void reset()
  {
    setDefault();
  }

  void setTolerance(Float val, int which=-1)
  {
    if (which<0)
      {
	for (long j=0;j<ncon;j++) 
	  ubvec[j]=val;
      }
    else ubvec[which]=val;
  }

  void setType(PartitionType type_)
  {
    type=type_;
  }

  void setRefineType(RefinePartitionType type_)
  {
    refineType=type_;
  }
  
  ParmetisParams(MpiCommunication *mpiCom_=glb::mpiComWorld,int ncon_=1):
    ncon(ncon_),mpiCom(mpiCom_)
  {
    nParts=mpiCom->size();
    tpwgts=(Float*) calloc(ncon*nParts,sizeof(Float));
    ubvec=(Float*) calloc(ncon,sizeof(Float));
    setDefault();
  }

  ~ParmetisParams()
  {
    free(tpwgts);tpwgts=NULL;
    free(ubvec);ubvec=NULL;
    freeData();
  }

  void freeData()
  {
    if (owner)
      {
	free(xadj);xadj=NULL;
	free(adjncy);adjncy=NULL;
	free(vtxdist);vtxdist=NULL;
	free(xyz);xyz=NULL;
	free(vwgt);vwgt=NULL;
	free(adjwgt);adjwgt=NULL;
      }
    //setDefault();
  }

  void setDefault()
  {
    wgtflag=0; // No weight
    numflag=0; // C-style
    opt[0]=0;opt[1]=0;opt[2]=0; // options [use options,output timing,seed]
    //opt[3]= PARMETIS_PSR_COUPLED ; // for adaptiverepart
    opt[3]= 1; // = PARMETIS_PSR_COUPLED 
    //printf("%d %d %lg\n",nParts,ncon,1.0/double(nParts));
    for (int i=0;i<nParts*ncon;i++) 
      tpwgts[i]=Float(1.0)/nParts; // weight expected on each partition (1/nParts)
    for (long j=0;j<ncon;j++) 
      ubvec[j]=1.05; //tolerance on repartition
    //vsize=1;
    itr=1.0; // should be within [0.0001, 1000000.0], 1000.0 recommanded by parmetis as default
    // the higher itr, the less edge cuts
    owner=false;

    vwgt=NULL;
    adjwgt=NULL;
    xadj=NULL;
    adjncy=NULL;
    vtxdist=NULL;
    xyz=NULL;
    vsize=NULL;
    nDims=0;
  }

  bool refine(std::vector<Index> &partition)
  {
    unsigned long nNodes=vtxdist[mpiCom->rank()+1]-vtxdist[mpiCom->rank()];
    MPI_Comm com=mpiCom->getCom(); // communicator
    partition.resize(nNodes);
    //printf("Nnodes = %ld\n",nNodes);
    if (refineType == RefinePartitionTypeV::PH)
      {
	ParMETIS_V3_PartGeom(vtxdist,&nDims,xyz,&partition[0],&com);
	edgeCut=-1;
      }	    
    if (refineType == RefinePartitionTypeV::KWAY)
      {
	ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,
			     vwgt,adjwgt,&wgtflag,
			     &numflag,&ncon,&nParts,
			     tpwgts,ubvec,opt,
			     &edgeCut,&partition[0],&com);
      }
    if (refineType == RefinePartitionTypeV::REFINE_KWAY)
      {
	ParMETIS_V3_RefineKway(vtxdist,xadj,adjncy,
			       vwgt,adjwgt,&wgtflag,
			       &numflag,&ncon,&nParts,
			       tpwgts,ubvec,opt,
			       &edgeCut,&partition[0],&com);
      }
    else if (refineType == RefinePartitionTypeV::ADAPTIVE)
      {
	ParMETIS_V3_AdaptiveRepart(vtxdist,xadj,adjncy,
				   vwgt,adjwgt,vsize, &wgtflag,
				   &numflag,&ncon,&nParts,
				   tpwgts,ubvec,&itr,opt,
				   &edgeCut,&partition[0],&com);
      }
    return true;
  }

  bool partition(std::vector<Index> &partition)
  {
    unsigned long nNodes=vtxdist[mpiCom->rank()+1]-vtxdist[mpiCom->rank()];
    MPI_Comm com=mpiCom->getCom(); // communicator
    partition.resize(nNodes);
    if (type == PartitionTypeV::KWAY)
      {
	// Purely topological version
	// Slow but gives a good partition
	// Prefered version if you have enough memory

	ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,
			     vwgt,adjwgt,&wgtflag,
			     &numflag,&ncon,&nParts,
			     tpwgts,ubvec,
			     opt,&edgeCut,&partition[0],&com);

      }
    else if (type == PartitionTypeV::PH)
      {
	// Purely geometrical, using PH-curve
	// very fast and cheap 
	// Sometimes the partition is good, sometimes not ...

	ParMETIS_V3_PartGeom(vtxdist,&nDims,xyz,&partition[0],&com);
	edgeCut=-1;
      }	
    else if (type == PartitionTypeV::PH_KWAY)
      {
	// Hybrid method
	// Not sure if the quality is better or worse than KWAY
	// uses more memory, but faster
	   	
	ParMETIS_V3_PartGeomKway(vtxdist,xadj,adjncy,
				 vwgt,adjwgt,&wgtflag,
				 &numflag,&nDims,xyz,&ncon,&nParts,
				 tpwgts,ubvec,
				 opt,&edgeCut,&partition[0],&com);
      }
    else 
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Unknown partition type.\n");
	exit(-1);
      }

    return true;
  }

};

#include "../internal/namespace.footer"
#endif
