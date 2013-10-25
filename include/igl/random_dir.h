#ifndef IGL_RANDOM_DIR_H
#define IGL_RANDOM_DIR_H
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl
{
  // Generate a uniformly random unit direction in 3D, return as vector
  Eigen::Vector3d random_dir();
  // Generate n stratified uniformly random unit directions in 3d, return as rows
  // of an n by 3 matrix
  //
  // Inputs:
  //   n  number of directions
  // Return n by 3 matrix of random directions
  Eigen::MatrixXd random_dir_stratified(const int n);
}

#ifdef IGL_HEADER_ONLY
#  include "random_dir.cpp"
#endif

#endif
