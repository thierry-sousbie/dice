#ifndef __MEMORY_POOL_ITERATORS_HXX__
#define __MEMORY_POOL_ITERATORS_HXX__

#include "../iterableMemoryPoolPolicies.hxx"

#include "../../../internal/namespace.header"

namespace internal {
  template <class C, class V, class Policy>
  class MemoryPoolIteratorT;
}

#include "../../../internal/namespace.footer"

// specialisations
#include "memoryPoolIteratorAlternate.hxx"
#include "memoryPoolIteratorBlocks.hxx"
#include "memoryPoolIteratorSortedBlocks.hxx"

#endif
