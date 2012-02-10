#ifndef IGL_ALL_PAIRS_DISTANCES_H
#define IGL_ALL_PAIRS_DISTANCES_H
#include "igl_inline.h"

namespace igl
{
  // ALL_PAIRS_DISTANCES compute distances between each point i in V and point j
  // in U
  // 
  // D = all_pairs_distances(V,U)
  // 
  // Templates:
  //   Mat  matrix class like MatrixXd
  // Inputs:
  //   V  #V by dim list of points
  //   U  #U by dim list of points
  //   squared  whether to return squared distances
  // Outputs:
  //   D  #V by #U matrix of distances, where D(i,j) gives the distance or
  //     squareed distance between V(i,:) and U(j,:)
  // 
  template <typename Mat>
  IGL_INLINE void all_pairs_distances(
    const Mat & V,
    const Mat & U,
    const bool squared, 
    Mat & D);
}

#ifdef IGL_HEADER_ONLY
#  include "all_pairs_distances.cpp"
#endif

#endif
