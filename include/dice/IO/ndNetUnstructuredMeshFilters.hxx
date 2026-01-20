#ifndef __DICE_NDNET_UNSTRUCTURED_MESH_DEFAULT_FILTERS_HXX__
#define __DICE_NDNET_UNSTRUCTURED_MESH_DEFAULT_FILTERS_HXX__

#include <functional>

/**
 * @file
 * @brief  Definition of a few default filters for NDnet unstructured mesh outputs
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup IO
 *   \{
 */

namespace IO
{

  /**
   * \class NDnetFilter_dataThresholdT
   * \brief A filter for NDnetUnstructuredMeshWriterT that filter cells by thresholding
   * with respect to a given cell data
   * \tparam M The mesh class
   * \tparam CellType The type of cells to filter (M::Vertex or M::Simplex)
   * \tparam Compare  Binary function that accepts two double elements as arguments : the
   * first is the vakue associated to a cell and the second a threshold. The function
   * returns true if the cell passed the test, false otherwise.
   */
  template <class M, class CellType, class Compare = std::less<double>>
  class NDnetFilter_dataThresholdT
  {
  public:
    typedef CellType Cell;

    /** \brief constructor.
     *  \param mesh A pointer to the mesh
     *  \param name A string representing the data to compare (must a valid cell functor)
     *  \param threshold The filter threshold
     *  \param cmp The comparison functor (see Compare)
     */
    NDnetFilter_dataThresholdT(const M *mesh, const std::string &name,
                               double threshold, Compare cmp = Compare())
    {
      functor = mesh->getSimplexFunctorPtr(name);
      th = threshold;
      compare = cmp;
    }
    /*
    NDnetFilter_dataThresholdT(const M *mesh, const std::string &name,
            double threshold)
    {
      simplexFunctor=mesh->getSimplexFunctorPtr(cellData);
      th=threshold;
      compare=Compare();
    }
    */
    Cell *operator()(Cell *cell) const
    {
      if (compare(functor(cell), th))
        return cell;
      else
        return static_cast<Cell *>(NULL);
    }

  private:
    typename M::SimplexFunctor *functor;
    double th;
    Compare compare;
  };

  template <class M, class Compare>
  class NDnetFilter_dataThresholdT<M, typename M::Vertex, Compare>
  {
  public:
    typedef typename M::Vertex Cell;

    NDnetFilter_dataThresholdT(const M *mesh, const std::string &name,
                               double threshold, Compare cmp = Compare())
    {
      functor = mesh->getVertexFunctorPtr(name);
      th = threshold;
      compare = cmp;
    }
    /*
    NDnetFilter_dataThresholdT(const M *mesh, const std::string &name,
             double threshold)
    {
      vertexFunctor=mesh->getVertexFunctorPtr();
      th=threshold;
      compare=Compare();
    }
    */

    Cell *operator()(Cell *cell) const
    {
      if (compare(functor(cell), th))
        return cell;
      else
        return static_cast<Cell *>(NULL);
    }

  private:
    typename M::VertexFunctor *functor;
    double th;
    Compare compare;
  };

} // namespace IO

/** \}*/
#include "../internal/namespace.footer"
#endif
