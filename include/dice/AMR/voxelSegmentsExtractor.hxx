#ifndef __VOXEL_SEGMENTS_EXTRACTOR_HXX__
#define __VOXEL_SEGMENTS_EXTRACTOR_HXX__

#include "../dice_globals.hxx"

/**
 * @file 
 * @brief  A class to extract edges from a set of voxels
 * @author Thierry Sousbie
 */
#include "../internal/namespace.header"

/** \addtogroup AMR
 *   \{
 */

/**
 * \class LocalAmrGridProjectorT
 *
 */

template <int ND, typename C=double>
class voxelSegmentsExtractorT
{
public:
  typedef C Coord;
  
  static int NDIM = ND;
  

  voxelSegmentsExtractorT()
  {}

  ~voxelSegmentsExtractorT()
  {}

  void reset()
  {

  }


private:
  typedef std::list<Coord> CoordList;

  std::vector< CoordList > coordListPool;

  std::list< CoordList* > coord[ND];
  std::list< SegList > coord[ND];
  

}

/** \}*/
#include "../internal/namespace.footer"
#endif
