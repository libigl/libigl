#ifndef IGL_REORDER_H
#define IGL_REORDER_H
#include "igl_inline.h"
#include <vector>
// For size_t
#include <stddef.h>
#include <cstdlib>

namespace igl
{
  // Act like matlab's Y = X[I] for std vectors
  // where I contains a vector of indices so that after,
  // Y[j] = X[I[j]] for index j
  // this implies that Y.size() == I.size()
  // X and Y are allowed to be the same reference
  template< class T >
  IGL_INLINE void reorder(
    const std::vector<T> & unordered,
    std::vector<size_t> const & index_map,
    std::vector<T> & ordered);
}

#ifdef IGL_HEADER_ONLY
#  include "reorder.cpp"
#endif

#endif
