#include "polar_svd.h"
#include <Eigen/SVD>

// From Olga's CGAL mentee's ARAP code
template<typename Mat>
IGL_INLINE void igl::polar_svd(const Mat& A, Mat& R, Mat& T)
{
  typedef Eigen::Matrix<typename Mat::Scalar,Mat::RowsAtCompileTime,1> Vec;
  Eigen::JacobiSVD<Mat> svd;
  svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
  const Mat& u = svd.matrixU();
  const Mat& v = svd.matrixV();
  const Vec& w = svd.singularValues();
  R = u*v.transpose();
  T = v*w.asDiagonal()*v.adjoint();
}

IGL_INLINE void igl::polar_svd(const Eigen::Matrix3f& A, Eigen::Matrix3f& R, Eigen::Matrix3f& T)
{
  typedef Eigen::Matrix<Eigen::Matrix3f::Scalar,3,1> Vec;
  Eigen::JacobiSVD<Eigen::Matrix3f> svd;
  svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
  const Eigen::Matrix3f& u = svd.matrixU();
  const Eigen::Matrix3f& v = svd.matrixV();
  const Vec& w = svd.singularValues();
  R = u*v.transpose();
  T = v*w.asDiagonal()*v.adjoint();
}

// Clang is giving an annoying warning inside Eigen
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wconstant-logical-operand"
#endif
IGL_INLINE void igl::polar_svd(const Eigen::Matrix2f& A, Eigen::Matrix2f& R, Eigen::Matrix2f& T)
{
  typedef Eigen::Matrix<Eigen::Matrix2f::Scalar,2,1> Vec;
  Eigen::JacobiSVD<Eigen::Matrix2f> svd;
  svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
  const Eigen::Matrix2f& u = svd.matrixU();
  const Eigen::Matrix2f& v = svd.matrixV();
  const Vec& w = svd.singularValues();
  R = u*v.transpose();
  T = v*w.asDiagonal()*v.adjoint();
}
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::polar_svd<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
#endif
