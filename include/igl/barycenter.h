#ifndef IGL_BARYCENTER_H
#define IGL_BARYCENTER_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // BARYCENTER
  //
  // B = barycenter(V,F)
  //
  // Compute the barycenter of every triangle
  //
  // Inputs:
  //   V #V x dim matrix of vertex coordinates
  //   F #F x simplex_size  matrix of indices of triangle corners
  // Output:
  //   BC a #F x dim matrix of 3d vertices
  IGL_INLINE void barycenter(
      const Eigen::MatrixXd & V,
      const Eigen::MatrixXi & F,
      Eigen::MatrixXd & BC);
}

#ifdef IGL_HEADER_ONLY
#  include "barycenter.cpp"
#endif

#endif
