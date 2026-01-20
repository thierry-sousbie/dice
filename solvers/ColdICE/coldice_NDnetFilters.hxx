#ifndef __COLDICE_NDNET_FILTERS_HXX__
#define __COLDICE_NDNET_FILTERS_HXX__

#include <functional>
#include <string.h>

namespace internal
{
  template <class M>
  class NDnetFilter_domainIndexedBaseT
  {
  public:
    static const int NDIM = M::NDIM;
    static const int fac = (NDIM == 2) ? 2 : 6;

    typedef typename M::Simplex Simplex;
    NDnetFilter_domainIndexedBaseT(M *mesh, int *initMeshResolution,
                                   const char *domainIndex = "domainIndex")
    {
      std::copy_n(initMeshResolution, NDIM, resolution);
      functor = mesh->getSimplexFunctorPtr(domainIndex);
      pad[0] = 1;
      for (int i = 0; i < NDIM; ++i)
        pad[i + 1] = pad[i] * resolution[i];
    }
    /*
      Simplex *operator()(Simplex *s)
      {
      int coords[NDIM];
      index2coords((*functor)(s),coords);
      }
    */
  protected:
    const typename M::SimplexFunctor *functor;
    int resolution[NDIM];
    long pad[NDIM + 1];

    void index2coords(const Simplex *s, int *coords) const
    {
      index2coords((long)(*functor)(s), coords);
    }

    void index2coords(long index, int *coords) const
    {
      index /= fac;
      coords[NDIM - 1] = index / pad[NDIM - 1];
      for (int i = NDIM - 2; i >= 0; --i)
      {
        index -= coords[i + 1] * pad[i + 1];
        coords[i] = index / pad[i];
      }
    }

    long shrunkIndex(const Simplex *s, const int factor[NDIM]) const
    {
      return shrunkIndex((long)(*functor)(s), factor);
    }

    long shrunkIndex(long index, const int factor[NDIM]) const
    {
      int coords[NDIM];
      index2coords(index, coords);
      for (int i = 0; i < NDIM; ++i)
        coords[i] /= factor[i];
      int newIndex = coords[0];
      int delta = 1;
      for (int i = 1; i < NDIM; ++i)
      {
        delta *= (resolution[i - 1] / factor[i - 1]);
        newIndex += coords[i] * delta;
      }
      return newIndex;
    }
  };
}

template <class M, class Compose = std::logical_and<bool>>
class NDnetFilter_subsetsT : public internal::NDnetFilter_domainIndexedBaseT<M>
{
public:
  typedef internal::NDnetFilter_domainIndexedBaseT<M> Base;

  static const int NDIM = Base::NDIM;
  typedef typename Base::Simplex Simplex;

  NDnetFilter_subsetsT(M *mesh, int *initMeshResolution,
                       const char *patterns[NDIM],
                       const char *domainIndex = "domainIndex") : Base(mesh, initMeshResolution, domainIndex)
  {
    initPattern(patterns);
  }

  NDnetFilter_subsetsT(M *mesh, int *initMeshResolution,
                       const char *domainIndex = "domainIndex") : Base(mesh, initMeshResolution, domainIndex)
  {
    const char *patterns[] = {"*", "_*", "*_"};
    initPattern(patterns);
  }

  Simplex *operator()(Simplex *s) const
  {
    int coords[NDIM];
    Base::index2coords(s, coords);
    bool result = compose(pat_[0][coords[0]], pat_[1][coords[1]]);
    for (int i = 2; i < NDIM; ++i)
      result = compose(result, pat_[i][coords[i]]);
    // printf("id=%ld: %d(%d) %d(%d) %d(%d) => %d\n",
    // 	   (long)(*Base::functor)(s),
    // 	   coords[0],(int)pat_[0][coords[0]],
    // 	   coords[1],(int)pat_[1][coords[1]],
    // 	   coords[2],(int)pat_[2][coords[2]],
    // 	   (int)result);

    if (result)
      return s;
    else
      return static_cast<Simplex *>(NULL);
  }

private:
  void initPattern(const char *patterns[NDIM])
  {
    for (int dim = 0; dim < NDIM; ++dim)
    {
      int ci = 0;
      int len = strlen(patterns[dim]);
      pat_[dim].reserve(Base::resolution[dim]);
      for (int i = 0; i < Base::resolution[dim]; ++i, ++ci)
      {
        if (ci >= len)
          ci = 0;
        pat_[dim].push_back((patterns[dim][ci] - '_') != 0);
      }
    }
  }
  Compose compose;
  std::vector<bool> pat_[NDIM];
};

// 2D
template <class M, int ND = M::NDIM>
class NDnetFilter_lagrangianLinesT : public internal::NDnetFilter_domainIndexedBaseT<M>
{
public:
  typedef internal::NDnetFilter_domainIndexedBaseT<M> Base;

  static const int NDIM = Base::NDIM;
  typedef typename Base::Simplex Simplex;
  typedef typename M::FacetHandle Handle;

  NDnetFilter_lagrangianLinesT(M *mesh, int *initMeshResolution,
                               const int delta[NDIM],
                               const char *domainIndex = "domainIndex") : Base(mesh, initMeshResolution, domainIndex)
  {
    std::copy_n(delta, NDIM, delta_);
  }

  Handle operator()(const Handle &h) const
  {
    if (Base::shrunkIndex(h->getSimplex(), delta_) !=
        Base::shrunkIndex(h->getOppositeSimplex(), delta_))
      return h;
    else
      return Handle();
  }

private:
  int delta_[NDIM];
};

// 3D
template <class M>
class NDnetFilter_lagrangianLinesT<M, 3> : public internal::NDnetFilter_domainIndexedBaseT<M>
{
public:
  typedef internal::NDnetFilter_domainIndexedBaseT<M> Base;

  static const int NDIM = Base::NDIM;
  typedef typename Base::Simplex Simplex;
  typedef typename M::SegmentHandle Handle;

  NDnetFilter_lagrangianLinesT(M *mesh, int *initMeshResolution,
                               const int delta[NDIM],
                               const char *domainIndex = "domainIndex") : Base(mesh, initMeshResolution, domainIndex)
  {
    std::copy_n(delta, NDIM, delta_);
  }

  Handle operator()(const Handle &h) const
  {
    auto it = h->getCirculator();
    const auto it_end = it;

    int val[3];
    int nVal = 1;

    val[0] = (int)(Base::shrunkIndex(*it, delta_));
    for (++it; (it != it_end); ++it)
    {
      int id = (int)(Base::shrunkIndex(*it, delta_));
      if ((std::distance(val, std::find(val, val + nVal, id)) == nVal) && (nVal < 3))
        val[nVal++] = id;
    }
    if (nVal > 2)
      return h;
    else
      return Handle();
  }

private:
  int delta_[NDIM];
};

#endif
