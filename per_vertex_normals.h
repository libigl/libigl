#ifndef IGL_PER_VERTEX_NORMALS_H
#define IGL_PER_VERTEX_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>
// Note: So for this only computes normals per vertex as uniformly weighted
// averages of incident triangle normals. It would be nice to support more or
// all of the methods here:
// "A comparison of algorithms for vertex normal computation"
namespace igl
{
  // Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  // Output:
  //   N  #V by 3 eigen Matrix of mesh vertex 3D normals
  IGL_INLINE void per_vertex_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & N);
}

#ifdef IGL_HEADER_ONLY
#  include "per_vertex_normals.cpp"
#endif

#endif
