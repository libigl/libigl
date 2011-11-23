#ifndef IGL_ADJACENCY_MATRIX_H
#define IGL_ADJACENCY_MATRIX_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  // Constructs the graph adjacency matrix  of a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F  #F by dim list of mesh faces (must be triangles)
  // Outputs: 
  //   A  max(F) by max(F) cotangent matrix, each row i corresponding to V(i,:)
  //
  // Example:
  //   // Mesh in (V,F)
  //   Eigen::SparseMatrix<double> A;
  //   adjacency_matrix(F,A);
  //   // sum each row 
  //   SparseVector<double> Asum;
  //   sum(A,1,Asum);
  //   // Convert row sums into diagonal of sparse matrix
  //   SparseMatrix<double> Adiag;
  //   diag(Asum,Adiag);
  //   // Build uniform laplacian
  //   SparseMatrix<double> U;
  //   U = A-Adiag;
  //
  // See also: edges, cotmatrix, diag
  template <typename T>
  inline void adjacency_matrix(
    const Eigen::MatrixXi & F, 
    Eigen::SparseMatrix<T>& A);
}

// Implementation
#include "verbose.h"

template <typename T>
inline void igl::adjacency_matrix(
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

#endif
