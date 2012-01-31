#include "print_ijv.h"

#include "find.h"
#include <iostream>

template <typename T>
IGL_INLINE void igl::print_ijv(
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

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
