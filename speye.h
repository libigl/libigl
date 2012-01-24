#ifndef IGL_SPEYE_H
#define IGL_SPEYE_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace igl
{
  // Builds an m by n sparse identity matrix like matlab's speye function
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   m  number of rows
  //   n  number of cols
  // Outputs:
  //   I  m by n sparse matrix with 1's on the main diagonal
  template <typename T>
  inline void speye(const int n,const int m, Eigen::SparseMatrix<T> & I);
  // Builds an n by n sparse identity matrix like matlab's speye function
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   n  number of rows and cols
  // Outputs:
  //   I  n by n sparse matrix with 1's on the main diagonal
  template <typename T>
  inline void speye(const int n, Eigen::SparseMatrix<T> & I);
}

// Implementation

template <typename T>
inline void igl::speye(const int m, const int n, Eigen::SparseMatrix<T> & I)
{
  // size of diagonal
  int d = (m<n?m:n);
  I = Eigen::SparseMatrix<T>(m,n);
  I.reserve(d);
  for(int i = 0;i<d;i++)
  {
    I.insert(i,i) = 1.0;
  }
  I.finalize();
}

template <typename T>
inline void igl::speye(const int n, Eigen::SparseMatrix<T> & I)
{
  return igl::speye(n,n,I);
}

#endif
