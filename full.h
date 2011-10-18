#ifndef IGL_FULL_H
#define IGL_FULL_H
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // Convert a sparsematrix into a full one
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Input:
  //   A  m by n sparse matrix
  // Output:
  //   B  m by n dense/full matrix
  template <typename T>
  inline void full(
    const Eigen::SparseMatrix<T> & A,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B);
}

// Implementation

template <typename T>
inline void igl::full(
  const Eigen::SparseMatrix<T> & A,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
  B = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(A.rows(),A.cols());
  // Iterate over outside
  for(int k=0; k<A.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (A,k); it; ++it)
    {
      B(it.row(),it.col()) = it.value();
    }
  }
}
#endif
