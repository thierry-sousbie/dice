#ifndef __COLDICE_CELL_DATA_HXX__
#define __COLDICE_CELL_DATA_HXX__

#include "config.h"

#ifndef D_DONT_USE_TRACERS

#define TRACER_REFINE_METHOD \
  dice::simplexRefineDataPolicy::RegressionWithTracerAndSegTracers_tracer
#define SEGTRACERS_REFINE_METHOD \
  dice::simplexRefineDataPolicy::RegressionWithTracerAndSegTracers_segTracers
#define VERTEX_REFINE_COORDS_METHOD \
  dice::vertexRefineCoordsPolicy::QuadraticRegressionWithTracers
#define SPLIT_LONGEST_EDGE_DEFAULT 0

#else

#define TRACER_REFINE_METHOD dice::simplexRefineDataPolicy::Barycenter
#define SEGTRACERS_REFINE_METHOD dice::simplexRefineDataPolicy::SegmentsBarycenter
#define VERTEX_REFINE_COORDS_METHOD dice::vertexRefineCoordsPolicy::MidPoint
#define SPLIT_LONGEST_EDGE_DEFAULT 1

#endif

#ifndef NO_SIMPLEX_TRACERS
#define SIMPLEX_TRACERS_INDEX 1
#else
#define SIMPLEX_TRACERS_INDEX 0
#endif

#if D_PER_SIMPLEX_INVARIANT
#define PER_SIMPLEX_INVARIANT_INDEX (SIMPLEX_TRACERS_INDEX +1)
#else
#define PER_SIMPLEX_INVARIANT_INDEX (SIMPLEX_TRACERS_INDEX)
#endif



#include <dice/tools/helpers/helpers_macros.hxx>
#include <dice/mesh/cellData/cellDataMacros.hxx>
#include <dice/mesh/cellData/cellDataFunctors_interface.hxx>

template <int D> class ColdiceVertexData;
template <int D> class ColdiceSimplexData;

template <>
class ColdiceVertexData<2>
{
public: 
  INIT_CELL_DATA(ColdiceVertexData<2>,Vertex);

  typedef dice::VertexDataElementT<0,2,double,
				   dice::vertexInitDataPolicy::Coords,
				   dice::vertexRefineDataPolicy::MidPoint> InitCoords;
  typedef dice::VertexDataElementT<1,1,double,
				   dice::vertexInitDataPolicy::Copy,
				   dice::vertexRefineDataPolicy::Average> ProjectedDensity;
  // typedef dice::VertexDataElementT<2,1,double,
  // 				   dice::vertexInitDataPolicy::Copy,
  // 				   dice::vertexRefineDataPolicy::Average> DisplacementField;
  InitCoords initCoords;  
  ProjectedDensity projectedDensity;
  // DisplacementField displacementField;
};
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<2>,InitCoords,initCoords)
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<2>,ProjectedDensity,projectedDensity)
// DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<2>,DisplacementField,displacementField)

template <>
class ColdiceVertexData<4>
{
public:
  INIT_CELL_DATA(ColdiceVertexData<4>,Vertex);
  
  typedef dice::VertexDataElementT<0,4,double,
				   dice::vertexInitDataPolicy::Coords,
				   dice::vertexRefineDataPolicy::MidPoint> InitCoords;
  typedef dice::VertexDataElementT<1,1,double,
				   dice::vertexInitDataPolicy::Copy,
				   dice::vertexRefineDataPolicy::Average> ProjectedDensity;
  // typedef dice::VertexDataElementT<2,2,double,
  // 				   dice::vertexInitDataPolicy::Copy,
  // 				   dice::vertexRefineDataPolicy::Average> DisplacementField;
  InitCoords initCoords; 
  ProjectedDensity projectedDensity;
  // DisplacementField displacementField;
};
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<4>,InitCoords,initCoords)
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<4>,ProjectedDensity,projectedDensity)
// DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<4>,DisplacementField,displacementField)

template <>
class ColdiceVertexData<6>
{
public:
  INIT_CELL_DATA(ColdiceVertexData<6>,Vertex);

  typedef dice::VertexDataElementT<0,6,double,
				   dice::vertexInitDataPolicy::Coords,
				   dice::vertexRefineDataPolicy::MidPoint> InitCoords;
  typedef dice::VertexDataElementT<1,1,double,
				   dice::vertexInitDataPolicy::Copy,
				   dice::vertexRefineDataPolicy::Average> ProjectedDensity;
  // typedef dice::VertexDataElementT<2,3,double,
  // 				   dice::vertexInitDataPolicy::Copy,
  // 				   dice::vertexRefineDataPolicy::Average> DisplacementField;
  InitCoords initCoords;  
  ProjectedDensity projectedDensity;
  // DisplacementField displacementField;
};
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<6>,InitCoords,initCoords)
DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<6>,ProjectedDensity,projectedDensity)
// DECLARE_CELL_DATA_ELEMENT( ColdiceVertexData<6>,DisplacementField,displacementField)

template <>
class ColdiceSimplexData<4>
{
public:
  INIT_CELL_DATA(ColdiceSimplexData<4>,Simplex);

  typedef dice::SimplexDataElementT<0,1,double,
				    dice::simplexInitDataPolicy::ProjectedVolumeWeighted,
				    dice::simplexRefineDataPolicy::HalfHalf,
				    dice::simplexCoarsenDataPolicy::Add> 
  Mass;
  
  typedef dice::SimplexDataElementT<1,2,double,
				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::Copy,
				    dice::simplexCoarsenDataPolicy::Average> 
  ProjectedDensityGradient;
  
  // The coarsen policy should be changed ... 
  typedef dice::SimplexDataElementT<2,4*3,double,
				    dice::simplexInitDataPolicy::SegTracers,
				    SEGTRACERS_REFINE_METHOD,
				    dice::simplexCoarsenDataPolicy::Average,
				    D_DUMP_TRACERS>
  SegTracers;

  // The coarsen policy should be changed ... 
#ifndef NO_SIMPLEX_TRACERS
  typedef dice::SimplexDataElementT<3+SIMPLEX_TRACERS_INDEX,4,double,
   				    dice::simplexInitDataPolicy::Barycenter,
   				    TRACER_REFINE_METHOD,
   				    dice::simplexCoarsenDataPolicy::Average> 
  Tracer;
#endif

  
#if D_PER_SIMPLEX_INVARIANT
  typedef dice::SimplexDataElementT<3+PER_SIMPLEX_INVARIANT_INDEX,1,double,		    
   				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::HalfHalf,
				    dice::simplexCoarsenDataPolicy::Add> 
  InvariantThreshold;
#endif
  /*
  typedef dice::SimplexDataElementT<4,1,long,
				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::Copy,
				    dice::simplexCoarsenDataPolicy::Keep>  
  DomainIndex;
  */
  Mass mass; 
  ProjectedDensityGradient projectedDensityGradient;
  SegTracers segTracers;
#ifndef NO_SIMPLEX_TRACERS
  Tracer tracer;
#endif
#if D_PER_SIMPLEX_INVARIANT
  InvariantThreshold invariantThreshold;
#endif
  //DomainIndex domainIndex;
};
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<4>,Mass,mass)
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<4>,ProjectedDensityGradient,projectedDensityGradient)
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<4>,SegTracers,segTracers)
#ifndef NO_SIMPLEX_TRACERS
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<4>,Tracer,tracer)
#endif
//DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<4>,DomainIndex,domainIndex)

template <>
class ColdiceSimplexData<6>
{
public:
  INIT_CELL_DATA(ColdiceSimplexData<6>,Simplex);

  typedef dice::SimplexDataElementT<0,1,double,
				    dice::simplexInitDataPolicy::ProjectedVolumeWeighted,
				    dice::simplexRefineDataPolicy::HalfHalf,
				    dice::simplexCoarsenDataPolicy::Add> 
  Mass;
  
  typedef dice::SimplexDataElementT<1,3,double,
				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::Copy,
				    dice::simplexCoarsenDataPolicy::Average> 
  ProjectedDensityGradient;

  // The coarsen policy should be changed ... 
  typedef dice::SimplexDataElementT<2,6*6,double,
				    dice::simplexInitDataPolicy::SegTracers,
				    SEGTRACERS_REFINE_METHOD,
				    dice::simplexCoarsenDataPolicy::Average,
				    D_DUMP_TRACERS>  
  SegTracers;

  // The coarsen policy should be changed ... 
#ifndef NO_SIMPLEX_TRACERS
  typedef dice::SimplexDataElementT<3+SIMPLEX_TRACERS_INDEX,6,double,
   				    dice::simplexInitDataPolicy::Barycenter,
   				    TRACER_REFINE_METHOD,
   				    dice::simplexCoarsenDataPolicy::Average> 
  Tracer;
#endif

#if D_PER_SIMPLEX_INVARIANT
  typedef dice::SimplexDataElementT<3+PER_SIMPLEX_INVARIANT_INDEX,1,double,		    
   				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::HalfHalf,
				    dice::simplexCoarsenDataPolicy::Add> 
  InvariantThreshold;
#endif
  
  /*
  typedef dice::SimplexDataElementT<4,1,long,
				    dice::simplexInitDataPolicy::Copy,
				    dice::simplexRefineDataPolicy::Copy,
				    dice::simplexCoarsenDataPolicy::Keep>  
  DomainIndex;
  */
  Mass mass; 
  ProjectedDensityGradient projectedDensityGradient;
  SegTracers segTracers;
#ifndef NO_SIMPLEX_TRACERS
  Tracer tracer;
#endif
#if D_PER_SIMPLEX_INVARIANT
  InvariantThreshold invariantThreshold;
#endif
  // DomainIndex domainIndex;
};
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<6>,Mass,mass)
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<6>,ProjectedDensityGradient,projectedDensityGradient)
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<6>,SegTracers,segTracers)
#ifndef NO_SIMPLEX_TRACERS
DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<6>,Tracer,tracer)
#endif
//DECLARE_CELL_DATA_ELEMENT(ColdiceSimplexData<6>,DomainIndex,domainIndex)

#endif
