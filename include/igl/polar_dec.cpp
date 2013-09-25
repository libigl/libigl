#include "polar_dec.h"
#include <Eigen/Eigenvalues>

#include "polar_svd.h"
#ifdef _WIN32
#else
#  include <fenv.h>
#endif
// You will need the development version of Eigen which is > 3.0.3
// You can determine if you have computeDirect by issuing
//   grep -r computeDirect path/to/eigen/*
#define EIGEN_HAS_COMPUTE_DIRECT

// From Olga's CGAL mentee's ARAP code
template<typename Mat>
IGL_INLINE void igl::polar_dec(const Mat& A, Mat& R, Mat& T)
{
#ifdef EIGEN_HAS_COMPUTE_DIRECT
 typedef typename Mat::Scalar Scalar;
 typedef Eigen::Matrix<typename Mat::Scalar,3,1> Vec;

 const Scalar th = std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());

 Eigen::SelfAdjointEigenSolver<Mat> eig;
 feclearexcept(FE_UNDERFLOW);
 eig.computeDirect(A.transpose()*A);
 if(fetestexcept(FE_UNDERFLOW) || eig.eigenvalues()(0)/eig.eigenvalues()(2)<th)
   return polar_svd(A,R,T);

 Vec S = eig.eigenvalues().cwiseSqrt();

 T = eig.eigenvectors() * S.asDiagonal() * eig.eigenvectors().transpose();
 R = A  * eig.eigenvectors() * S.asDiagonal().inverse()
        * eig.eigenvectors().transpose();

 if(std::abs(R.squaredNorm()-3.) > th)
   return polar_svd(A,R,T);
#else
  return polar_svd(A,R,T);
#endif
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::polar_dec<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
#endif
