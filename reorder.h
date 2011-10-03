#ifndef IGL_REORDER_H
#define IGL_REORDER_H

namespace igl
{
  // Act like matlab's Y = X[I] for std vectors
  // where I contains a vector of indices so that after,
  // Y[j] = X[I[j]] for index j
  // this implies that Y.size() == I.size()
  // X and Y are allowed to be the same reference
  template< class T >
  inline void reorder(
    const std::vector<T> & unordered,
    std::vector<size_t> const & index_map,
    std::vector<T> & ordered);
}

// Implementation

// This implementation is O(n), but also uses O(n) extra memory
template< class T >
inline void igl::reorder(
  const std::vector<T> & unordered,
  std::vector<size_t> const & index_map,
  std::vector<T> & ordered)
{
  // copy for the reorder according to index_map, because unsorted may also be
  // sorted
  std::vector<T> copy = unordered;
  ordered.resize(index_map.size());
  for(int i = 0; i<(int)index_map.size();i++)
  {
    ordered[i] = copy[index_map[i]];
  }
}
#endif

