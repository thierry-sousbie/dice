#ifndef __LOCAL_AMR_GRID_THREADED_ALLOCATOR_HXX__
#define __LOCAL_AMR_GRID_THREADED_ALLOCATOR_HXX__

#include "../tools/memory/memoryPool.hxx"

#include "../internal/namespace.header"

namespace localAmrGrid
{

  template <class V, int NT>
  class ThreadedAllocatorT
  {
  protected:
    typedef V Voxel;
    static const int N_THREADS = NT;

    template <class V>
    struct VoxelGroupT
    {
      static const long NDIM = V::NDIM;
      static const long LEVEL = L;
      static const long N_CHILD = (1 << NDIM);

      Voxel data[N_CHILD_TOTAL];
    };

  public:
    typedef VoxelGroupT<V> VoxelGroup;

    void pop(VoxelGroup **vg)
    {
      voxelGroupPool[omp_get_thread_num()].pop(vg);
    }

    void pop(VoxelGroup **vg, int threadIndex)
    {
      voxelGroupPool[threadIndex].pop(vg);
    }

  protected:
    typedef MemoryPoolT<VoxelGroupT<Voxel>> Pool;
    Pool pool[NT];

  private:
    /*
    namespace internal {

      template <int N, int COUNT>
      struct PoolTypeList
      {
  static const int LEVEL = N;
  typedef MemoryPoolT< VoxelGroup<D,N> > Pool;
  typedef PoolTypeList<N+1,COUNT> Next;

  Pool pool;
      };

      template <int COUNT>
      struct PoolTypeList<COUNT,COUNT>
      {
  static const int LEVEL = COUNT;
  typedef MemoryPoolT< VoxelGroup<D,COUNT> > Pool;
  typedef PoolTypeList<COUNT,COUNT> Next;

  Pool pool;
      };

    }
    */
  };

}

/** \}*/
#include "../internal/namespace.footer"
#endif
