#ifndef __MESH_ITERATORS_MACROS_HXX__
#define __MESH_ITERATORS_MACROS_HXX__

#include "../tools/helpers/helpers_macros.hxx"

#define FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads) \
  _Pragma(STRINGIFY( CONCAT_ONE_ARGUMENT(omp parallel for schedule(dynamic,1) num_threads,nThreads) ))

#define FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads) \
  _Pragma(STRINGIFY( CONCAT_ONE_ARGUMENT(omp parallel for schedule(static,1) num_threads,nThreads) ))


// Batch macros
#define FOREACH_BATCH_SIMPLEX(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);	\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_iterator itVar=mesh->simplexBegin(iVar,foreach_loop_count);\
    const simplexPtr_iterator itVar##_end=mesh->simplexEnd();\

#define FOREACH_BATCH_SIMPLEX_LG(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_LG_iterator itVar=mesh->simplexLGBegin(iVar,foreach_loop_count);\
    const simplexPtr_LG_iterator itVar##_end=mesh->simplexLGEnd();\

#define FOREACH_BATCH_SIMPLEX_LGS(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_LGS_iterator itVar=mesh->simplexLGSBegin(iVar,foreach_loop_count);\
    const simplexPtr_LGS_iterator itVar##_end=mesh->simplexLGSEnd();\

#define FOREACH_BATCH_VERTEX(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_iterator itVar=mesh->vertexBegin(iVar,foreach_loop_count);\
    const vertexPtr_iterator itVar##_end=mesh->vertexEnd();\

#define FOREACH_BATCH_VERTEX_LG(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_LG_iterator itVar=mesh->vertexLGBegin(iVar,foreach_loop_count);\
    const vertexPtr_LG_iterator itVar##_end=mesh->vertexLGEnd();\

#define FOREACH_BATCH_VERTEX_LGS(mesh,nThreads,nBatches,iVar,itVar)\
  {const int foreach_loop_count = (nThreads==1)?1:(nThreads*nBatches);\
  FOREACH_BATCH_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_LGS_iterator itVar=mesh->vertexLGSBegin(iVar,foreach_loop_count);\
    const vertexPtr_LGS_iterator itVar##_end=mesh->vertexLGSEnd();\

#define FOREACH_BATCH_END }}

// thread macros (i.e. =1 batch)
#define FOREACH_THREAD_SIMPLEX(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_iterator itVar=mesh->simplexBegin(iVar,nThreads);\
    const simplexPtr_iterator itVar##_end=mesh->simplexEnd();\

#define FOREACH_THREAD_SIMPLEX_LG(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_LG_iterator itVar=mesh->simplexLGBegin(iVar,nThreads);\
    const simplexPtr_LG_iterator itVar##_end=mesh->simplexLGEnd();\

#define FOREACH_THREAD_SIMPLEX_LGS(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    simplexPtr_LGS_iterator itVar=mesh->simplexLGSBegin(iVar,nThreads);\
    const simplexPtr_LGS_iterator itVar##_end=mesh->simplexLGSEnd();\

#define FOREACH_THREAD_VERTEX(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_iterator itVar=mesh->vertexBegin(iVar,nThreads);\
    const vertexPtr_iterator itVar##_end=mesh->vertexEnd();\

#define FOREACH_THREAD_VERTEX_LG(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_LG_iterator itVar=mesh->vertexLGBegin(iVar,nThreads);\
    const vertexPtr_LG_iterator itVar##_end=mesh->vertexLGEnd();\

#define FOREACH_THREAD_VERTEX_LGS(mesh,nThreads,iVar,itVar)\
  {const int foreach_loop_count = nThreads;\
  FOREACH_THREAD_MESH_CELL_OPENMP_ARG(nThreads)\
  for (int iVar=0;iVar<foreach_loop_count;++iVar){\
    vertexPtr_LGS_iterator itVar=mesh->vertexLGSBegin(iVar,nThreads);\
    const vertexPtr_LGS_iterator itVar##_end=mesh->vertexLGSEnd();\

#define FOREACH_THREAD_END }}

#endif


