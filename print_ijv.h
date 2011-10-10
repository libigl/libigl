#ifndef IGL_PRINT_IJV_H
#define IGL_PRINT_IJV_H
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // Prints a 3 column matrix representing [I,J,V] = find(X). That is, each
  // row is the row index, column index and value for each non zero entry. Each
  // row is printed on a new line
  //
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Input:
  //   X  m by n matrix whose entries are to be sorted
  //   offset  optional offset for I and J indices {0}
  template <typename T>
  inline void print_ijv(
    const Eigen::SparseMatrix<T>& X, 
    const int offset=0);
}

// Implementation
#include "find.h"
#include <iostream>

template <typename T>
inline void igl::print_ijv(
  const Eigen::SparseMatrix<T>& X,
  const int offset)
{
  Eigen::Matrix<int,Eigen::Dynamic,1> I;
  Eigen::Matrix<int,Eigen::Dynamic,1> J;
  Eigen::Matrix<T,Eigen::Dynamic,1> V;
  igl::find(X,I,J,V);
  // Concatenate I,J,V
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> IJV(I.size(),3);
  IJV.col(0) = I.cast<T>();
  IJV.col(1) = J.cast<T>();
  IJV.col(2) = V;
  // Offset
  if(offset != 0)
  {
    IJV.col(0).array() += offset;
    IJV.col(1).array() += offset;
  }
  std::cout<<IJV;
}

#endif
