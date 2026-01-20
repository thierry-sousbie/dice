#ifndef __LOCAL_AMR_GRID_VOXEL_HXX__
#define __LOCAL_AMR_GRID_VOXEL_HXX__

#include "./mpiStructs/mpiExchangeStruct_voxelData.hxx"

/**
 * @file
 * @brief Voxel structure desgined to be used by LocalAmrGridT
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup AMR
 *   \{
 */

// #define LOCAL_AMR_GRID_VOXEL_STORE_PARENTS

/**
 * \class LocalAmrGridVoxelT
 * \brief Voxel structure designed to be used by LocalAmrGridT
 *
 * \tparam G  AMR grid type (most probably LocalAmrGridT)
 */
template <class G>
class LocalAmrGridVoxelT
{

public:
  typedef LocalAmrGridVoxelT<G> MyType;
  typedef G Grid;

  typedef typename G::ICoord ICoord;
  typedef typename G::Data Data;

  static const int NDIM = G::NDIM;
  static const int NVERT = (1 << G::NDIM);
  static const int NSEG = ((1 << G::NDIM) * G::NDIM) / 2;

  static const int NNEI = (1 << G::NDIM);
  static const int NSEGNEI = (1 << (G::NDIM - 1));
  static const int NVERTNEI = (1 << G::NDIM);

  typedef mpiExchangeStruct::Mpi_VoxelDataT<MyType> MpiExchangeData;

  LocalAmrGridVoxelT(char level_ = 0)
  {
    level = level_;
    child = NULL;
    data = Data();
    childId = -1;
  }

  //! reset the content of a voxel, but keeps its index
  void setEmpty()
  {
    child = NULL;
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
    parent = NULL;
#endif
    data = Data();
    childId = -1;
    /*
    level=0;
    flags=0;
    childId=-1;
    */
  }

  char getLevel() const
  {
    return level;
  }

  ICoord getIndex() const
  {
    return index;
  }

  MyType *getChild(int quadrant) const
  {
    return &child[quadrant];
  }

  bool isLeaf() const
  {
    return (child == NULL);
  }

  bool isRoot() const
  {
    return (level == 0);
  }

  /** return in which quadrant (octant in 3D) of this voxel
   * a voxel with index \a otherIndex falls
   */
  int getQuadrant(ICoord otherIndex) const
  {
    // return G::getQuadrantFromRef(index,otherIndex);
    return G::getQuadrantAtLevel(otherIndex, level);
  }

  int getQuadrant(MyType *other) const
  {
    return getQuadrant(other->index);
  }

  MyType *child;
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
  MyType *parent;
#endif
  ICoord index;
  Data data;

  unsigned int flags;
  char level;
  char childId;
  unsigned short poolId;

  void setIndex(ICoord index_)
  {
    index = index_;
  }

  MyType *refineCritical(Grid *grid)
  {
    /*
    MyType childBuffer[G::CHILDREN_COUNT];
    for (int i=0;i<G::CHILDREN_COUNT;++i)
      {
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
  tmp[i].parent  = this;
#endif
  childBuffer[i].level   = level+1;
  childBuffer[i].childId = i;
  childBuffer[i].index   = G::computeChildIndex(index,level,i);
      }

#pragma omp critical
    {
      if (isLeaf())
  {
    MyType *tmp;
    grid->popVoxelGroup(&tmp,level+1);
    memcpy(tmp,childBuffer,sizeof(MyType)*G::CHILDREN_COUNT);
    // Child == NULL indicates a leaf, so we do that last !
    child=tmp;
  }
    }
*/

#pragma omp critical
    {
      if (isLeaf())
      {
        MyType *tmp;
        grid->popVoxelGroup(&tmp, level + 1);
        for (int i = 0; i < G::CHILDREN_COUNT; ++i)
        {
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
          tmp[i].parent = this;
#endif
          tmp[i].level = level + 1;
          tmp[i].childId = i;
          tmp[i].index = G::computeChildIndex(index, level, i);
        }
        child = tmp;
      }
    }

    return this;
  }

  MyType *refine(Grid *grid)
  {
    if (isLeaf())
    {
      grid->popVoxelGroup(&child, level + 1);
      for (int i = 0; i < G::CHILDREN_COUNT; ++i)
      {
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
        child[i].parent = this;
#endif
        child[i].level = level + 1;
        child[i].childId = i;
        child[i].index = G::computeChildIndex(index, level, i);
      }
    }
    return this;
  }

  MyType *refine(Grid *grid, int level)
  {
    if (level > getLevel())
    {
      refine(grid);
      if (level > getLevel() + 1)
        for (int i = 0; i < G::CHILDREN_COUNT; ++i)
          child[i].refine_rec(grid, level);
    }
    return this;
  }

  void resetFlags()
  {
    flags = 0;
  }

  void resetVerticesFlags()
  {
    flags &= ~((1 << NVERT) - 1);
  }

  void setVertexFlag(int i)
  {
    flags |= (1 << i);
  }

  /** \brief returns whether the flag associated with vertex \a i is set. This usually
   *  means that this vertex belongs to the voxel.
   *  \param i index of the vertex to test (0 <= i < (1<<NDIM))
   *  \return true if the flag is set, false otherwise.
   */
  int getVertexFlag(int i) const
  {
    return flags & (1 << i);
  }

  /** \brief returns the (1<<NDIM) flags associated with the vertices.
   *  if this is 0, then the voxel does not own any of its vertices
   *  \return true if the flag is set, false otherwise.
   */
  int getVertexFlags() const
  {
    return flags & ((1 << NVERT) - 1);
  }

  void resetSegmentsFlags()
  {
    flags &= ((1 << NVERT) - 1);
  }

  void setSegmentFlag(int vertexId, int dim)
  {
    vertexId &= ~(1 << dim);
    flags |= (1 << dim) << (vertexId * G::NDIM + NVERT);
  }

  int getSegmentFlags() const
  {
    return flags >> NVERT;
  }

  /*
  int getSegmentFlag(int vertexId)
  {
    return (flags>>(vertexId*G::NDIM+NVERT))&((1<<G::NDIM)-1);
  }

  int getSegmentFlag(int vertexId, int dim)
  {
    return (flags>>(vertexId*G::NDIM+NVERT))&(1<<dim);
  }
  */
  void print(const G *grid, const std::string &txt)
  {
    double center[G::NDIM];
    grid->index2CenterCoords(index, center);
    if (G::NDIM == 2)
    {
      printf("%s Voxel %ld (%d child %d parents, flags=%ud) at level %d centered on (%f,%f) : %g\n",
             txt.c_str(),
             index, 8 * (child != NULL),
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
             parent != NULL,
#else
             -1 /*parent!=NULL*/,
#endif
             flags,
             level, center[0], center[1], data);
    }
    else if (G::NDIM == 3)
    {
      printf("%s Voxel %ld (%d child %d parents, flags=%ud) at level %d centered on (%f,%f,%f) : %g\n",
             txt.c_str(),
             index, 8 * (child != NULL),
#ifdef LOCAL_AMR_GRID_VOXEL_STORE_PARENTS
             parent != NULL,
#else
             -1 /*parent!=NULL*/,
#endif
             flags,
             level, center[0], center[1], center[2], data);
    }
  }

private:
  void refine_rec(Grid *grid, int level)
  {
    refine(grid);
    if (level > getLevel() + 1)
      for (int i = 0; i < G::CHILDREN_COUNT; ++i)
        child[i].refine_rec(grid, level);
  }

  /*
  void print(const G *grid, const std::string &txt)
  {
    double center[3];
    grid->index2CenterCoords(index,center);
    printf("%s Voxel %ld (%d child %d parents, flags=%ud) at level %d centered on (%f,%f,%f) : %g\n",
     txt.c_str(),
     index,8*(child!=NULL),parent!=NULL,flags,
     level,center[0],center[1],center[2],data);
  }
  */
};

/** \}*/
#include "../internal/namespace.footer"
#endif
