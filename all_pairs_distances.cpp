#include "all_pairs_distances.h"

template <typename Mat>
IGL_INLINE void igl::all_pairs_distances(
  const Mat & V,
  const Mat & U,
  const bool squared,
  Mat & D)
{
  // dimension should be the same
  assert(V.cols() == U.cols());
  // resize output
  D.resize(V.rows(),U.rows());
  for(int i = 0;i<V.rows();i++)
  {
    for(int j=0;j<U.rows();j++)
    {
      D(i,j) = (V.row(i)-U.row(j)).array().pow(2).sum();
      if(!squared)
      {
        D(i,j) = sqrt(D(i,j));
      }
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
