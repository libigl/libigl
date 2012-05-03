#include "adjacency_matrix.h"

#include "verbose.h"

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

template <typename T>
IGL_INLINE void igl::adjacency_matrix(
  const Eigen::MatrixXi & F, 
  Eigen::SparseMatrix<T>& A)
{
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> 
    dyn_A(F.maxCoeff()+1, F.maxCoeff()+1);
  dyn_A.reserve(6*(F.maxCoeff()+1));

  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,(j+1)%F.cols());
      dyn_A.coeffRef(s, d) = 1;
      dyn_A.coeffRef(d, s) = 1;
    }
  }

  A = Eigen::SparseMatrix<T>(dyn_A);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
