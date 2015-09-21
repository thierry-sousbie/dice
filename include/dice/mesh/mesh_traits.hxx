#ifndef __MESH_TRAITS_HXX__
#define __MESH_TRAITS_HXX__

#include <limits>

#include "../mesh/simplexType.hxx"
#include "../mesh/cellData/cellDataEmpty.hxx"
#include "../mesh/basicVertexCoordsPolicies.hxx"

#include "../tools/types/globalIdentity.hxx"

#include "../geometry/boundaryType.hxx"

/**
 * @file 
 * @brief  Traits of meshT class.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup MESH
 *   \{
 */

/** 
 * \class MeshTraitsT 
 * \brief Traits of meshT class.
 * \tparam D the number of dimensions of the mesh
 * \tparam DW the number of dimensions of the embedding space
 * \tparam BT boundary conditions (see BoundaryType)
 * \tparam OnRefineCoordsPolicy policy defining how to computed the coordinates of a new vertex when a simplex is refined. See vertexRefineCoordsPolicy.
 * \tparam VD Data type to store with each vertex
 * \tparam SD Data type to store with each simplex
 */

template <int D, int DW, int BT = BoundaryType::NONE,
	  template <class,class,class,class,class> class OnRefineCoordsPolicy = 
	  vertexRefineCoordsPolicy::MidPoint,	  
	  class VD = VertexDataEmpty,
	  class SD = SimplexDataEmpty> 
class MeshTraitsT
{
public:

  static const int NDIM = D;
  static const int NDIM_W = DW;
  static const int MAX_MPI_SIZE = 10000;    
  static const int SIMPLEX_TYPE = simplexType::VerticesOnly;  
  static const int BOUNDARY_TYPE = BT;
  static const int WORLD_BOUNDARY_TYPE = BoundaryType::NONE;  
  
  typedef double Coord;
  
  // Must be able to index the number of simplices on a single node
#ifdef USELONGINT
  typedef unsigned long LocalIndex; 
#else
  typedef unsigned int LocalIndex; 
#endif
  typedef unsigned long long GlobalIndex; 
  typedef GlobalIdentityT<GlobalIndex,MAX_MPI_SIZE> GlobalIdentity;

  typedef unsigned char SimplexFlag;
  typedef unsigned char SegmentFlag;
  typedef unsigned char VertexFlag;  

  typedef VD VertexData;
  typedef SD SimplexData;
  
  static const GlobalIndex GLOBAL_INDEX_INVALID;
  static const GlobalIndex GLOBAL_INDEX_MAX;
  static const LocalIndex  LOCAL_INDEX_INVALID;
  static const LocalIndex  LOCAL_INDEX_MAX;

  template <class M, class SEG, class V, class S>
  static long getSimplexRefineBufferSize()
  {
    return OnRefineCoordsPolicy<Coord,M,SEG,V,S>::SIMPLEX_BUFFER_SIZE;    
  }

  template <class M, class SEG, class V, class S>
  static void refineCoords(const M *mesh, const SEG &seg, V* newVertex, 
			   S * const *simplices, int nSimplices, void *buffer) 
  {
    OnRefineCoordsPolicy<Coord,M,SEG,V,S>::refine(newVertex->getCoordsPtr(),
						  mesh,seg,newVertex,
						  simplices,nSimplices,buffer);
  }
};

template <int D, int DW, int BT, template <class,class,class,class,class> class ORCP, class VD,class SD> 
const typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GlobalIndex
MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GLOBAL_INDEX_INVALID = 
  std::numeric_limits<typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GlobalIndex>::max();

template <int D, int DW, int BT, template <class,class,class,class,class> class ORCP,class VD,class SD> 
const typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GlobalIndex
MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GLOBAL_INDEX_MAX = 
  std::numeric_limits<typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::GlobalIndex>::max()-1;

template <int D, int DW, int BT, template <class,class,class,class,class> class ORCP,class VD,class SD> 
const typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LocalIndex
MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LOCAL_INDEX_INVALID = 
  std::numeric_limits<typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LocalIndex>::max();

template <int D, int DW, int BT, template <class,class,class,class,class> class ORCP,class VD,class SD> 
const typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LocalIndex
MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LOCAL_INDEX_MAX = 
  std::numeric_limits<typename MeshTraitsT<D,DW,BT,ORCP,VD,SD>::LocalIndex>::max()-1;

/** \}*/
#include "../internal/namespace.footer"
#endif
