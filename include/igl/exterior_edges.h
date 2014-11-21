#ifndef IGL_EXTERIOR_EDGES_H
#define IGL_EXTERIOR_EDGES_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // EXTERIOR_EDGES Determines boundary "edges" and also edges with an
  // odd number of occurances where seeing edge (i,j) counts as +1 and seeing
  // the opposite edge (j,i) counts as -1
  //
  // Inputs:
  //   F  #F by simplex_size list of "faces"
  // Outputs:
  //   E  #E by simplex_size-1  list of exterior edges
  //
  IGL_INLINE void exterior_edges(
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & E);
  // Inline version
  IGL_INLINE Eigen::MatrixXi exterior_edges( const Eigen::MatrixXi & F);
}
#ifndef IGL_STATIC_LIBRARY
#  include "exterior_edges.cpp"
#endif

#endif
