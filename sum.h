#ifndef IGL_SUM_H
#define IGL_SUM_H
#include <Eigen/Sparse>

namespace igl
{
  // Ideally, this becomes a super overloaded function that works with sparse
  // and dense matrices like the matlab sum function

  // Sum the columns or rows of a sparse matrix
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   X  m by n sparse matrix
  //   dim  dimension along which to sum (1 or 2)
  // Output:
  //   S  n-long sparse vector (if dim == 1) 
  //   or
  //   S  m-long sparse vector (if dim == 2)
  template <typename T>
  inline void sum(
    const Eigen::SparseMatrix<T>& X, 
    const int dim,
    Eigen::SparseVector<T>& S);
}

template <typename T>
inline void igl::sum(
  const Eigen::SparseMatrix<T>& X, 
  const int dim,
  Eigen::SparseVector<T>& S)
{
  // dim must be 2 or 1
  assert(dim == 1 || dim == 2);
  // Get size of input
  int m = X.rows();
  int n = X.cols();
  // resize output
  if(dim==1)
  {
    S = Eigen::SparseVector<T>(n);
  }else
  {
    S = Eigen::SparseVector<T>(m);
  }

  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      if(dim == 1)
      {
        S.coeffRef(it.col()) += it.value();
      }else
      {
        S.coeffRef(it.row()) += it.value();
      }
    }
  }

}
#endif
