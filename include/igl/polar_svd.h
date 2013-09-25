#ifndef IGL_POLAR_SVD
#define IGL_POLAR_SVD
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // Computes the polar decomposition (R,T) of a matrix A using SVD singular value decomposition
  // Inputs:
  //   A  3 by 3 matrix to be decomposed
  // Outputs:
  //   R  3 by 3 rotation matrix part of decomposition
  //   T  3 by 3 stretch matrix part of decomposition
  //
  IGL_INLINE void polar_svd(const Eigen::Matrix3f& A, Eigen::Matrix3f& R, Eigen::Matrix3f& T);
  IGL_INLINE void polar_svd(const Eigen::Matrix2f& A, Eigen::Matrix2f& R, Eigen::Matrix2f& T);
  template<typename Mat>
  IGL_INLINE void polar_svd(const Mat& A, Mat& R, Mat& T);
}
#ifdef IGL_HEADER_ONLY
#  include "polar_svd.cpp"
#endif
#endif
