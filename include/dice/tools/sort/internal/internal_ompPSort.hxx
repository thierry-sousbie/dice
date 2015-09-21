#ifndef __INTERNAL_OPEN_MP_PARALLEL_SORT__
#define __INTERNAL_OPEN_MP_PARALLEL_SORT__

#include <algorithm> 

#include "../../../internal/namespace.header"

namespace internal {

  template< class T, class Compare>
  void ompPSort(T* data, unsigned long len, int maxThreadLen, const Compare &comp)
  {
    if (len <= maxThreadLen)
      std::sort(data, data + len, comp);
    else
      {
	unsigned long half_len = len/2;
#pragma omp task shared(comp) untied
	ompPSort<T,Compare>(data, half_len, maxThreadLen, comp);
	  
#pragma omp task shared(comp) untied
	ompPSort<T,Compare>(data + half_len, len - half_len, maxThreadLen, comp); 

#pragma omp taskwait
	std::inplace_merge(data, data + half_len, data + len, comp);
      }
  }

}

#include "../../../internal/namespace.footer"
#endif
