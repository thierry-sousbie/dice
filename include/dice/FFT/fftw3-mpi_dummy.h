#ifndef __FFTW3_MPI_DUMMY_H__
#define __FFTW3_MPI_DUMMY_H__

#ifdef HAVE_FFTW3_MPI
#include <fftw3-mpi.h>

bool fftwNeedsPadding(void *src, void *dst)
{
  return true;
}

template <class T>
void fftwTransposedSwap(T &a, T&b)
{
  std::swap(a,b);
}

#else

#ifndef FFTW_MPI_TRANSPOSED_OUT 
#define FFTW_MPI_TRANSPOSED_OUT 0
#endif

#ifndef FFTW_MPI_TRANSPOSED_IN
#define FFTW_MPI_TRANSPOSED_IN 0 
#endif

bool fftwNeedsPadding(void *src, void *dst)
{
  return (src==dst);
}

template <class T>
void fftwTransposedSwap(T &a, T&b)
{}

ptrdiff_t fftw_mpi_local_size(int rnk, const ptrdiff_t *n, void* comm,
			      ptrdiff_t *local_n0, ptrdiff_t *local_0_start)
{
  ptrdiff_t result=1;
  *local_n0 = n[0];
  *local_0_start = 0;
  for (int i=0;i<rnk;++i)
    result*=n[i];
  //result*=(n[rnk-1]/2+1);

  return result;
}

ptrdiff_t fftw_mpi_local_size_transposed(int rnk, const ptrdiff_t *n, void* comm,
					 ptrdiff_t *local_n0, ptrdiff_t *local_0_start,
					 ptrdiff_t *local_n1, ptrdiff_t *local_1_start)
{
  ptrdiff_t result=1;
  *local_n0 = n[0];
  *local_0_start = 0;
  *local_n1 = n[1];
  *local_1_start = 0;

  for (int i=0;i<rnk;++i)
    result*=n[i];
  //result*=(n[rnk-1]/2+1);

  return result;
}

template <typename R, typename C>
fftw_plan fftw_mpi_plan_dft_r2c(int rnk, const ptrdiff_t *n, R *in, C *out, 
				void *comm, unsigned flags)
{
  std::vector<int> nn(n,n+rnk);
  return fftw_plan_dft_r2c(rnk,&nn[0],in,out,flags);
}

template <typename R, typename C>
fftw_plan fftw_mpi_plan_dft_c2r(int rnk, const ptrdiff_t *n, C *in, R *out, 
				void *comm, unsigned flags)
{
  std::vector<int> nn(n,n+rnk);
  return fftw_plan_dft_c2r(rnk,&nn[0],in,out,flags);
}


#endif
#endif
