#ifndef IGL_MEDIAN_H
#define IGL_MEDIAN_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Compute the median of an eigen vector
  //
  // Inputs:
  //   V  #V list of unsorted values
  // Outputs:
  //   m  median of those values
  // Returns true on success, false on failure
  IGL_INLINE bool median(const Eigen::VectorXd & V, double & m);
}

#ifdef IGL_HEADER_ONLY
#  include "median.cpp"
#endif

#endif
