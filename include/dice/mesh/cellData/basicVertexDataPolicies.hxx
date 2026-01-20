#ifndef __BASIC_VERTEX_DATA_POLICIES_HXX__
#define __BASIC_VERTEX_DATA_POLICIES_HXX__

#include "../../internal/namespace.header"

namespace vertexInitDataPolicy
{

  template <typename T, int N, class M, class V, class DT>
  struct Zero
  {
    static void init(T *result, const M *mesh, const V *vertex, const DT *initVal)
    {
      for (int i = 0; i < N; ++i)
        result[i] = T();
    }
  };

  template <typename T, int N, class M, class V, class DT>
  struct Copy
  {
    static void init(T *result, const M *mesh, const V *vertex, const DT *initVal)
    {
      for (int i = 0; i < N; ++i)
        result[i] = (initVal == NULL) ? T() : initVal[i];
    }
  };

  template <typename T, int N, class M, class V, class DT>
  struct Coords
  {
    static void init(T *result, const M *mesh, const V *vertex, const DT *initVal)
    {
      for (int i = 0; i < N; ++i)
        result[i] = vertex->getCoord(i);
    }
  };

}

namespace vertexRefineDataPolicy
{

  template <int I, typename T, int N, class M, class SEG, class V, class S>
  struct Dummy
  {
    static void refine(T *result, const M *mesh, const SEG &seg, const V *newVertex,
                       S *const *simplices, int nSimplices) {}
  };

  template <int I, typename T, int N, class M, class SEG, class V, class S>
  struct MidPoint
  {
    static void refine(T *result, const M *mesh, const SEG &seg, const V *newVertex,
                       S *const *simplices, int nSimplices)
    {
      T *v0 = seg->getVertex(0)->template getDataElementPtr<I>()->getPointer();
      T *v1 = seg->getVertex(1)->template getDataElementPtr<I>()->getPointer();
      mesh->getGeometry()->midPointCoords(v0, v1, result);
    }
  };

  template <int I, typename T, int N, class M, class SEG, class V, class S>
  struct Average
  {
    static void refine(T *result, const M *mesh, const SEG &seg, const V *newVertex,
                       S *const *simplices, int nSimplices)
    {
      T *v0 = seg->getVertex(0)->template getDataElementPtr<I>()->getPointer();
      T *v1 = seg->getVertex(1)->template getDataElementPtr<I>()->getPointer();
      for (int i = 0; i < N; ++i)
        result[i] = (v0[i] + v1[i]) * 0.5;
    }
  };

}

#include "../../internal/namespace.footer"
#endif
