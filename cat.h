#ifndef IGL_CAT_H
#define IGL_CAT_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace igl
{
  // If you're using Dense matrices you might be better off using the << operator

  // This is an attempt to act like matlab's cat function.

  // Perform concatenation of a two matrices along a single dimension
  // If dim == 1, then C = [A;B]. If dim == 2 then C = [A B]
  // 
  // Template:
  //   Scalar  scalar data type for sparse matrices like double or int
  //   Mat  matrix type for all matrices (e.g. MatrixXd, SparseMatrix)
  //   MatC  matrix type for ouput matrix (e.g. MatrixXd) needs to support
  //     resize
  // Inputs:
  //   A  first input matrix
  //   B  second input matrix
  //   dim  dimension along which to concatenate, 0 or 1
  // Outputs:
  //   C  output matrix
  //   
  template <typename Scalar>
  inline void cat(
      const int dim, 
      const Eigen::SparseMatrix<Scalar> & A, 
      const Eigen::SparseMatrix<Scalar> & B, 
      Eigen::SparseMatrix<Scalar> & C);
  template <typename Derived, class MatC>
  inline void cat(
    const int dim,
    const Eigen::MatrixBase<Derived> & A, 
    const Eigen::MatrixBase<Derived> & B,
    MatC & C);
  // Wrapper that returns C
  template <class Mat>
  inline Mat cat(const int dim, const Mat & A, const Mat & B);

  // Note: Maybe we can autogenerate a bunch of overloads D = cat(int,A,B,C),
  // E = cat(int,A,B,C,D), etc. 

  // Concatenate a "matrix" of blocks
  // C = [A0;A1;A2;...;An] where Ai = [A[i][0] A[i][1] ... A[i][m]];
  //
  // Inputs:
  //   A  a matrix (vector of row vectors)
  // Output:
  //   C
  template <class Mat>
  inline void cat(const std::vector<std::vector< Mat > > & A, Mat & C);
}

// Implementation

// Sparse matrices need to be handled carefully. Because C++ does not 
// Template:
//   Scalar  sparse matrix scalar type, e.g. double
template <typename Scalar>
inline void igl::cat(
    const int dim, 
    const Eigen::SparseMatrix<Scalar> & A, 
    const Eigen::SparseMatrix<Scalar> & B, 
    Eigen::SparseMatrix<Scalar> & C)
{
  assert(dim == 1 || dim == 2);
  using namespace Eigen;
  // Special case if B or A is empty
  if(A.size() == 0)
  {
    C = B;
    return;
  }
  if(B.size() == 0)
  {
    C = A;
    return;
  }

  DynamicSparseMatrix<Scalar, RowMajor> dyn_C;
  if(dim == 1)
  {
    assert(A.cols() == B.cols());
    dyn_C.resize(A.rows()+B.rows(),A.cols());
  }else if(dim == 2)
  {
    assert(A.rows() == B.rows());
    dyn_C.resize(A.rows(),A.cols()+B.cols());
  }else
  {
    fprintf(stderr,"cat.h: Error: Unsupported dimension %d\n",dim);
  }

  dyn_C.reserve(A.nonZeros()+B.nonZeros());

  // Iterate over outside of A
  for(int k=0; k<A.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename SparseMatrix<Scalar>::InnerIterator it (A,k); it; ++it)
    {
      dyn_C.coeffRef(it.row(),it.col()) += it.value();
    }
  }

  // Iterate over outside of B
  for(int k=0; k<B.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename SparseMatrix<Scalar>::InnerIterator it (B,k); it; ++it)
    {
      int r = (dim == 1 ? A.rows()+it.row() : it.row());
      int c = (dim == 2 ? A.cols()+it.col() : it.col());
      dyn_C.coeffRef(r,c) += it.value();
    }
  }

  C = SparseMatrix<Scalar>(dyn_C);
}

template <typename Derived, class MatC>
inline void igl::cat(
  const int dim,
  const Eigen::MatrixBase<Derived> & A, 
  const Eigen::MatrixBase<Derived> & B,
  MatC & C)
{
  assert(dim == 1 || dim == 2);
  // Special case if B or A is empty
  if(A.size() == 0)
  {
    C = B;
    return;
  }
  if(B.size() == 0)
  {
    C = A;
    return;
  }

  if(dim == 1)
  {
    assert(A.cols() == B.cols());
    C.resize(A.rows()+B.rows(),A.cols());
    C << A,B;
  }else if(dim == 2)
  {
    assert(A.rows() == B.rows());
    C.resize(A.rows(),A.cols()+B.cols());
    C << A,B;
  }else
  {
    fprintf(stderr,"cat.h: Error: Unsupported dimension %d\n",dim);
  }
}

template <class Mat>
inline Mat igl::cat(const int dim, const Mat & A, const Mat & B)
{
  assert(dim == 1 || dim == 2);
  Mat C;
  igl::cat(dim,A,B,C);
  return C;
}

template <class Mat>
inline void cat(const std::vector<std::vector< Mat > > & A, Mat & C)
{
  using namespace igl;
  using namespace std;
  // Start with empty matrix
  C.resize(0,0);
  for(typename vector<vector< Mat > >::const_iterator rit = A.begin(); rit != A.end(); rit++)
  {
    // Concatenate each row horizontally
    // Start with empty matrix
    Mat row(0,0);
    for(typename vector<vector< Mat > >::iterator cit = A.begin(); rit != A.end(); rit++)
    {
      row = cat(2,row,*cit);
    }
    // Concatenate rows vertically
    C = cat(1,C,row);
  }
}

#endif
