#ifndef IGL_REPMAT_H
#define IGL_REPMAT_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // Ideally this is a super overloaded function that behaves just like
  // matlab's repmat

  // Replicate and tile a matrix
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  m by n input matrix
  //   r  number of row-direction copies
  //   c  number of col-direction copies
  // Outputs:
  //   B  r*m by c*n output matrix
  //
  template <typename T,const int W, const int H>
  IGL_INLINE void repmat(
    const Eigen::Matrix<T,W,H> & A,
    const int r,
    const int c,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B);
  template <typename T>
  IGL_INLINE void repmat(
    const Eigen::SparseMatrix<T> & A,
    const int r,
    const int c,
    Eigen::SparseMatrix<T> & B);
}

#ifdef IGL_HEADER_ONLY
#  include "repmat.cpp"
#endif

#endif
