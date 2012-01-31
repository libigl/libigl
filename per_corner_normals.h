#ifndef IGL_PER_CORNER_NORMALS_H
#define IGL_PER_CORNER_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   corner_threshold  threshold in degrees on sharp angles
  // Output:
  //   CN  #F*3 by 3 eigen Matrix of mesh vertex 3D normals, where the normal
  //     for corner F(i,j) is at CN(i*3+j,:) 
  IGL_INLINE void per_corner_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const double corner_threshold,
    Eigen::MatrixXd & CN);
  // Other Inputs:
  //   FN  #F by 3 eigen Matrix of face normals
  //   VF  map from vertices to list of incident faces
  IGL_INLINE void per_corner_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const std::vector<std::vector<int> >& VF,
    const double corner_threshold,
    Eigen::MatrixXd & CN);
}

#ifdef IGL_HEADER_ONLY
#  include "per_corner_normals.cpp"
#endif

#endif
