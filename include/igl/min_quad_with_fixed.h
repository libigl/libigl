#ifndef IGL_MIN_QUAD_WITH_FIXED_H
#define IGL_MIN_QUAD_WITH_FIXED_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseExtra>

namespace igl
{
  template <typename T>
  struct min_quad_with_fixed_data;
  // MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C with
  // constraints that Z(known) = Y, optionally also subject to the constraints
  // Aeq*Z = Beq
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  n by n matrix of quadratic coefficients
  //   B  n by 1 column of linear coefficients
  //   known list of indices to known rows in Z
  //   Y  list of fixed values corresponding to known rows in Z
  //   Optional:
  //     Aeq  m by n list of linear equality constraint coefficients
  //     Beq  m by 1 list of linear equality constraint constant values
  //     pd flag specifying whether A(unknown,unknown) is positive definite
  // Outputs:
  //   data  factorization struct with all necessary information to solve
  //     using min_quad_with_fixed_solve
  // Returns true on success, false on error
  template <typename T>
  IGL_INLINE bool min_quad_with_fixed_precompute(
    const Eigen::SparseMatrix<T>& A,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & known,
    const Eigen::SparseMatrix<T>& Aeq,
    const bool pd,
    min_quad_with_fixed_data<T> & data
    );

  // Solves a system previously factored using min_quad_with_fixed_precompute
  // Inputs:
  //   data  factorization struct with all necessary precomputation to solve
  // Outputs:
  //   Z  n by cols solution
  // Returns true on success, false on error
  template <typename T>
  IGL_INLINE bool min_quad_with_fixed_solve(
    const min_quad_with_fixed_data<T> & data,
    const Eigen::Matrix<T,Eigen::Dynamic,1> & B,
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Y,
    const Eigen::Matrix<T,Eigen::Dynamic,1> & Beq,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Z);
}

template <typename T>
struct igl::min_quad_with_fixed_data
{
  // Size of original system: number of unknowns + number of knowns
  int n;
  // Whether A(unknown,unknown) is positive definite
  bool Auu_pd;
  // Whether A(unknown,unknown) is symmetric
  bool Auu_sym;
  // Indices of known variables
  Eigen::Matrix<int,Eigen::Dynamic,1> known;
  // Indices of unknown variables
  Eigen::Matrix<int,Eigen::Dynamic,1> unknown;
  // Indices of lagrange variables
  Eigen::Matrix<int,Eigen::Dynamic,1> lagrange;
  // Indices of unknown variable followed by Indices of lagrange variables
  Eigen::Matrix<int,Eigen::Dynamic,1> unknown_lagrange;
  // Matrix multiplied against Y when constructing right hand side
  Eigen::SparseMatrix<T> preY;
  // Tells whether system is sparse
  bool sparse;
  // Lower triangle of LU decomposition of final system matrix
  Eigen::SparseMatrix<T> L;
  // Upper triangle of LU decomposition of final system matrix
  Eigen::SparseMatrix<T> U;
  // Dense LU factorization
  Eigen::FullPivLU<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > lu;
};

#ifdef IGL_HEADER_ONLY
#  include "min_quad_with_fixed.cpp"
#endif

#endif
