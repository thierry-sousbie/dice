#ifndef __OPEN_MP_PARALLEL_SORT__
#define __OPEN_MP_PARALLEL_SORT__

#include "internal/internal_ompPSort.hxx"
/**
 * @file 
 * @brief A openMP parallel wrapper to std::sort
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 *   \{
 */

template< class T, class Compare>
void ompPSort(T *data, long len, int nThreads,const Compare &comp)
{
  if (nThreads<=1) std::sort(data,data+len,comp);
  else
    {
      long maxThreadLen = (len/nThreads)+1;
      if (maxThreadLen < 1024) maxThreadLen=1024;

#pragma omp parallel num_threads(nThreads)
      {
#pragma omp single
	internal::ompPSort<T,Compare>(data,len,maxThreadLen,comp);
      }
    }
}

template< class T > //, class Compare = std::less<T> >
void ompPSort(T *data, long len, int nThreads)
{
  std::less<T> cmp;
  ompPSort(data,len,nThreads,cmp);
}

template< class T, class Compare>
void ompPSort(T data_begin, T data_end, int nThreads,const Compare &comp)
{  
  unsigned long len = std::distance(data_begin,data_end);
  ompPSort<typename T::value_type,Compare>(&(*data_begin),len,nThreads,comp);
}

template< class T>
void ompPSort(T data_begin, T data_end, int nThreads)
{  
  std::less<typename T::value_type> cmp;
  ompPSort(data_begin,data_end,nThreads,cmp);
}


/** \}*/
#include "../../internal/namespace.footer"
#endif
