#ifndef IGL_SLICE_INTO_H
#define IGL_SLICE_INTO_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // Act like the matlab Y(row_indices,col_indices) = X
  // 
  // Inputs:
  //   X  xm by xn rhs matrix
  //   R  list of row indices
  //   C  list of column indices
  //   Y  ym by yn lhs matrix
  // Output:
  //   Y  ym by yn lhs matrix, same as input but Y(R,C) = X
  template <typename T>
  IGL_INLINE void slice_into(
    const Eigen::SparseMatrix<T>& X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::SparseMatrix<T>& Y);

  template <typename T, const int W, const int H>
  IGL_INLINE void slice_into(
    const Eigen::Matrix<T,W,H> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::Matrix<T,W,H> & Y);

  template <typename T>
  IGL_INLINE void slice_into(
    const Eigen::Matrix<T,Eigen::Dynamic,1> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    Eigen::Matrix<T,Eigen::Dynamic,1> & Y);
}

#ifdef IGL_HEADER_ONLY
#  include "slice_into.cpp"
#endif

#endif
