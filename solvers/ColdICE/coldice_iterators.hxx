#ifndef __COLDICE_ITERATORS_HXX__
#define __COLDICE_ITERATORS_HXX__

/**
 * @file
 * @brief  A local AMR grid designed to allow the projection of an unstructured mesh
 * @author Thierry Sousbie
 */

/** \addtogroup ColDICE
 *   \{
 */

#include <iterator>

/**
 * \class MeshAndTracersCoordsT
 * \brief A dummy container class used to iterate over all the vertices and
 * tracers coordinates defined in a mesh in one go. Each iterator dereferences to
 * a pointer to M::Coord.
 * \tparam M the mesh class
 */

template <class M>
class MeshAndTracersCoordsT
{
private:
  class m_iterator;

public:
  typedef M Mesh;
  typedef typename Mesh::Coord Coord;
  typedef m_iterator iterator;

  MeshAndTracersCoordsT(Mesh *m) : mesh(m) {}

  iterator begin(int delta = 0, int stride = 1)
  {
    return m_iterator(mesh->vertexLGSBegin(delta, stride),
                      mesh->vertexLGSEnd(),
                      mesh->simplexLGSBegin(delta, stride),
                      mesh->simplexLGSEnd());
  }

  iterator end(int delta = 0, int stride = 1)
  {
    return m_iterator(mesh->vertexLGSEnd(),
                      mesh->vertexLGSEnd(),
                      mesh->simplexLGSEnd(),
                      mesh->simplexLGSEnd());
  }

  iterator begin(int delta = 0, int stride = 1) const
  {
    return m_iterator(mesh->vertexLGSBegin(delta, stride),
                      mesh->vertexLGSEnd(),
                      mesh->simplexLGSBegin(delta, stride),
                      mesh->simplexLGSEnd());
  }

  iterator end(int delta = 0, int stride = 1) const
  {
    return m_iterator(mesh->vertexLGSEnd(),
                      mesh->vertexLGSEnd(),
                      mesh->simplexLGSEnd(),
                      mesh->simplexLGSEnd());
  }

private:
  Mesh *mesh;
  typedef typename Mesh::vertexPtr_LGS_iterator v_iterator;
  typedef typename Mesh::simplexPtr_LGS_iterator s_iterator;

  class m_iterator
  {
  private:
#ifndef NO_SIMPLEX_TRACERS
    static const int TRACERS_COUNT = Mesh::Simplex::NSEG + 1;
#else
    static const int TRACERS_COUNT = Mesh::Simplex::NSEG;
#endif

    static const int INDEX_TH = TRACERS_COUNT * Mesh::NDIM_W;

  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef m_iterator self_type;

    typedef typename Mesh::Coord *value_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef const value_type *const_pointer;
    typedef const value_type &const_reference;
    typedef long difference_type;

    m_iterator(v_iterator v_begin, v_iterator v_end,
               s_iterator s_begin, s_iterator s_end) : vit(v_begin), vEnd(v_end),
                                                       sit(s_begin), sEnd(s_end),
                                                       index(0)
    {
    }

    self_type &operator++()
    {
      if (vit != vEnd)
        ++vit;
      else
      {
        index += Mesh::NDIM_W;
        if (index == INDEX_TH)
        {
          index = 0;
          ++sit;
        }
      }

      return *this;
    }

    value_type operator*()
    {
      if (vit != vEnd)
        return (*vit)->getCoordsPtr();
#ifndef NO_SIMPLEX_TRACERS
      else if (index == (TRACERS_COUNT - 1) * Mesh::NDIM_W)
        return (*sit)->tracer.getPointer();
#endif
      else
        return &((*sit)->segTracers.getPointer()[index]);
    }

    value_type operator*() const
    {
      if (vit != vEnd)
        return (*vit)->getCoordsConstPtr();
#ifndef NO_SIMPLEX_TRACERS
      else if (index == (TRACERS_COUNT - 1) * Mesh::NDIM_W)
        return (*sit)->tracer.getConstPointer();
#endif
      else
        return &((*sit)->segTracers.getConstPointer()[index]);
    }

    bool operator==(const self_type &r) const
    {
      return (sit == r.sit) && (vit == r.vit) && (index == r.index);
    }

    bool operator!=(const self_type &r) const
    {
      return (sit != r.sit) || (vit != r.vit) || (index != r.index);
    }

  private:
    v_iterator vit;
    v_iterator vEnd;
    s_iterator sit;
    s_iterator sEnd;
    int index;
  };
};

/**
 * \class MeshAndTracersSliceCoordsT
 * \brief A dummy container class used to iterate over all the vertices and
 * tracers coordinates defined in a slice of a mesh in one go. The slice is defined as by
 * its origin and extent. Each iterator dereferences to a pointer to M::Coord.
 * \tparam M the mesh class
 */
template <class M>
class MeshAndTracersSliceCoordsT
{
public:
  typedef MeshAndTracersCoordsT<M> MeshAndTracersCoords;
  typedef M Mesh;
  typedef typename Mesh::Coord Coord;
  typedef typename std::vector<Coord *>::iterator iterator;
  typedef typename std::vector<Coord *>::const_iterator const_iterator;

  static const int NDIM = M::NDIM;

  template <class G>
  static int getSlicesCount(const G &grid, double maxAlloc)
  {
    // Permanently disabled for now as feature is not implemented yet !
    return 1;

    if (maxAlloc <= 0)
      return 1;
    double gridSize = sizeof(typename G::Data);
    for (int i = 0; i < G::NDIM; ++i)
      gridSize *= grid.getResolution(i);

    int nSlicesMax = gridSize / maxAlloc;
    if (nSlicesMax <= 1)
      nSlicesMax = 1;

    return nSlicesMax;
  }

  template <class G>
  MeshAndTracersSliceCoordsT(Mesh *m, const G &grid, double maxAlloc)
  {
    dim = 0;
    currentGroup = 0;
    // Find the first dimension along which the global grid is sliced
    // As we also want ot slice the mesh along that direction !
    for (int i = 0; i < NDIM; ++i)
    {
      if ((grid.getResolution(i) != grid.getLocalResolution(i)) && (dim < 0))
        dim = i;
      gridDim[i] = grid.getResolution(i);
    }

    x0 = grid.getVertexCoord(dim).front();
    delta = grid.getVertexCoord(dim).back() - x0;
    delta_inv = 1.0 / delta;

    MeshAndTracersCoords coordContainer(m);
    int minIndex = gridDim[dim];

    // Find the min and max index of the coordinates in the mesh
#pragma omp parallel for num_threads(dice::glb::num_omp_threads) reduction(min : minIndex)
    for (int j = 0; j < dice::glb::num_omp_threads; j++)
    {

      const auto it_end = coordContainer.end(j, dice::glb::num_omp_threads);
      for (auto it = coordContainer.begin(j, dice::glb::num_omp_threads); it != it_end; ++it)
      {
        int index = ((*it)[dim] - x0) * delta_inv;
        if (index < minIndex)
          minIndex = index;
      }
    }

    int nSlicesMax = getSlicesCount(grid, maxAlloc);
    groups.resize(nSlicesMax);
    // And compute coordinates group for fast and balanced openMP iteration
    int sliceExtent = (gridDim[dim] / nSlicesMax);
    if (sliceExtent * nSlicesMax < gridDim[dim])
      sliceExtent++;

#pragma omp parallel for num_threads(dice::glb::num_omp_threads)
    for (int j = 0; j < dice::glb::num_omp_threads; j++)
    {
      std::vector<std::vector<Coord *>> localGroups;
      localGroups.resize(nSlicesMax);

      const auto it_end = coordContainer.end(j, dice::glb::num_omp_threads);
      for (auto it = coordContainer.begin(j, dice::glb::num_omp_threads); it != it_end; ++it)
      {
        int index = ((*it)[dim] - x0) * delta_inv;
        int sliceIndex = (index - minIndex) / sliceExtent;
        localGroups[sliceIndex].push_back(*it);
      }

#pragma omp critical
      {
        for (int i = 0; i < groups.size(); ++i)
          groups[i].insert(groups[i].end(), localGroups[i].begin(), localGroups[i].end());
      }
    }

    // move empty groups to the end ...
    for (int i = 0, j = 0; j < groups.size(); ++j)
    {
      if (groups[j].size() != 0)
        groups[i++].swap(groups[j]);
    }
  }

  long getGroup() const
  {
    return currentGroup;
  }

  void setGroup(int group)
  {
    currentGroup = group;
  }

  long nGroups() const
  {
    return groups.size();
  }

  iterator begin(int delta = 0, int stride = 1)
  {
    long sz = groups[currentGroup].size() / stride;
    iterator it = groups[currentGroup].begin();
    std::advance(it, sz * delta);
    return it;
  }

  iterator end(int delta = 0, int stride = 1)
  {
    if (delta == (stride - 1))
      return groups[currentGroup].end();

    long sz = groups[currentGroup].size() / stride;
    iterator it = groups[currentGroup].begin();
    std::advance(it, sz * (delta + 1));
    return it;
  }

  const_iterator begin(int delta = 0, int stride = 1) const
  {
    long sz = groups[currentGroup].size() / stride;
    const_iterator it = groups[currentGroup].cbegin();
    std::advance(it, sz * delta);
    return it;
  }

  const_iterator end(int delta = 0, int stride = 1) const
  {
    if (delta == (stride - 1))
      return groups[currentGroup].cend();

    long sz = groups[currentGroup].size() / stride;
    const_iterator it = groups[currentGroup].cbegin();
    std::advance(it, sz * (delta + 1));
    return it;
  }

private:
  std::vector<std::vector<Coord *>> groups;

  int currentGroup;
  int dim;
  int nSlices;
  int gridDim[NDIM];
  double x0;
  double delta;
  double delta_inv;
};

/** \}*/
#endif
