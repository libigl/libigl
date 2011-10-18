#ifndef IGL_IS_SYMMETRIC_H
#define IGL_IS_SYMMETRIC_H
namespace igl
{
  // Returns true if the given matrix is symmetric
  // Inputs:
  //   A  m by m matrix
  // Returns true if the matrix is square and symmetric
  template <typename T>
  inline bool is_symmetric(const Eigen::SparseMatrix<T>& A);
}

// Implementation

template <typename T>
inline bool igl::is_symmetric(const Eigen::SparseMatrix<T>& A)
{
  if(A.rows() != A.cols())
  {
    return false;
  }
  Eigen::SparseMatrix<T> AT = A.transpose();
  Eigen::SparseMatrix<T> AmAT = A-AT;
  //// Eigen screws up something with LLT if you try to do
  //SparseMatrix<T> AmAT = A-A.transpose();
  //// Eigen crashes at runtime if you try to do
  // return (A-A.transpose()).nonZeros() == 0;
  return AmAT.nonZeros() == 0;
}
#endif
