#ifndef IGL_BOUNDING_BOX_DIAGONAL_H
#define IGL_BOUNDING_BOX_DIAGONAL_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Compute the length of the diagonal of a given meshes axis-aligned bounding
  // box
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  // Returns length of bounding box diagonal
  IGL_INLINE double bounding_box_diagonal(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F);
}

#ifdef IGL_HEADER_ONLY
#  include "bounding_box_diagonal.cpp"
#endif

#endif
