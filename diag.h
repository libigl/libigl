#ifndef IGL_DIAG_H
#define IGL_DIAG_H
#include <Eigen/Sparse>

namespace igl
{
  // Either extracts the main diagonal of a matrix as a vector OR converts a
  // vector into a matrix with vector along the main diagonal. Like matlab's
  // diag function

  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   X  an m by n sparse matrix
  // Outputs:
  //   V  a min(m,n) sparse vector
  template <typename T>
  void diag(
    const Eigen::SparseMatrix<T>& X, 
    Eigen::SparseVector<T>& V);
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //  V  a m sparse vector
  // Outputs:
  //  X  a m by m sparse matrix
  template <typename T>
  void diag(
    const Eigen::SparseVector<T>& V,
    Eigen::SparseMatrix<T>& X);
}

// Implementation
#include "verbose.h"

template <typename T>
void igl::diag(
  const Eigen::SparseMatrix<T>& X, 
  Eigen::SparseVector<T>& V)
{
  // Get size of input
  int m = X.rows();
  int n = X.cols();
  V = Eigen::SparseVector<T>((m>n?n:m));

  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      if(it.col() == it.row())
      {
        V.coeffRef(it.col()) += it.value();
      }
    }
  }
}

template <typename T>
void igl::diag(
  const Eigen::SparseVector<T>& V,
  Eigen::SparseMatrix<T>& X)
{
  // clear and resize output
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(V.size(),V.size());

  // loop over non-zeros
  for(typename SparseVector<T>::InnerIterator it(V); it; ++it)
  {
    dyn_X.coeffRef(it.index(),it.index()) += it.value();
  }

  X = Eigen::SparseMatrix<T>(dyn_X);
}
#endif
