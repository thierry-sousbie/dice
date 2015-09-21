#ifndef __PARMETIS_INTERFACE_HXX__
#define __PARMETIS_INTERFACE_HXX__

#ifdef USE_MPI
#include <parmetis.h>

#else

#include "../tools/MPI/myMpi.hxx"

typedef int idx_t;
typedef float real_t;

/*-------------------------------------------------------------------
 * API Introduced with Release 3.0 (current API) 
 *--------------------------------------------------------------------*/

int ParMETIS_V3_PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
			 idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, 
			 idx_t *ncon, idx_t *nparts, 
			 real_t *tpwgts, real_t *ubvec, idx_t *options, 
			 idx_t *edgecut, idx_t *part, 
			 MPI_Comm *comm){return 0;}

int  ParMETIS_V3_PartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
			      idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, 
			      idx_t *ndims, real_t *xyz, 
			      idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
			      idx_t *options, 
			      idx_t *edgecut, idx_t *part, MPI_Comm *comm){return 0;}

int  ParMETIS_V3_PartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, 
			  MPI_Comm *comm){return 0;}

int  ParMETIS_V3_RefineKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
			    idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, 
			    idx_t *ncon, idx_t *nparts, 
			    real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, 
			    idx_t *part, MPI_Comm *comm){return 0;}

int  ParMETIS_V3_AdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
				idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, 
				idx_t *numflag, idx_t *ncon, 
				idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
				real_t *ipc2redist, 
				idx_t *options, idx_t *edgecut, 
				idx_t *part, MPI_Comm *comm){return 0;}

int  ParMETIS_V3_Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, 
			   idx_t *ncommonnodes, idx_t **xadj, idx_t **adjncy, 
			   MPI_Comm *comm){return 0;}

int  ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
			      idx_t *wgtflag, idx_t *numflag, idx_t *ncon, 
			      idx_t *ncommonnodes, idx_t *nparts, 
			      real_t *tpwgts, real_t *ubvec, idx_t *options, 
			      idx_t *edgecut, idx_t *part, 
			      MPI_Comm *comm){return 0;}

int  ParMETIS_V3_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
			idx_t *options, idx_t *order, idx_t *sizes, 
			MPI_Comm *comm){return 0;}

int  ParMETIS_V32_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
			 idx_t *numflag, idx_t *mtype, idx_t *rtype, 
			 idx_t *p_nseps, idx_t *s_nseps,
			 real_t *ubfrac, idx_t *seed, idx_t *dbglvl, idx_t *order, 
			 idx_t *sizes, MPI_Comm *comm){return 0;}

int  ParMETIS_SerialNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
			   idx_t *options, idx_t *order, idx_t *sizes, 
			   MPI_Comm *comm){return 0;}


#endif

#endif
