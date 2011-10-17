#ifndef IGL_MIN_QUAD_WITH_FIXED_H
#define IGL_MIN_QUAD_WITH_FIXED_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
  bool min_quad_with_fixed_precompute(
    const Eigen::SparseMatrix<T>& A,
    const Eigen::MatrixXi & known,
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
  bool min_quad_with_fixed_solve(
    const min_quad_with_fixed_data<T> & data,
    const Eigen::Matrix<T,Eigen::Dynamic,1> & B,
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Y,
    const Eigen::Matrix<T,Eigen::Dynamic,1> & Beq,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Z);
}

// Implementation
#include <Eigen/SparseExtra>
#include <cassert>
#include <cstdio>
#include "slice.h"
#include "is_symmetric.h"

#include "find.h"
#include "sparse.h"
#include "lu_lagrange.h"

template <typename T>
struct igl::min_quad_with_fixed_data
{
  int n;
  bool Auu_pd;
  bool Auu_sym;
  Eigen::Matrix<int,Eigen::Dynamic,1> known;
  Eigen::Matrix<int,Eigen::Dynamic,1> unknown;
  Eigen::Matrix<int,Eigen::Dynamic,1> lagrange;
  Eigen::Matrix<int,Eigen::Dynamic,1> unknown_lagrange;
  SparseMatrix<T> preY;
  SparseMatrix<T> L;
  SparseMatrix<T> U;
};

template <typename T>
bool igl::min_quad_with_fixed_precompute(
  const Eigen::SparseMatrix<T>& A,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & known,
  const Eigen::SparseMatrix<T>& Aeq,
  const bool pd,
  igl::min_quad_with_fixed_data<T> & data
  )
{
  // number of rows
  int n = A.rows();
  // cache problem size
  data.n = n;

  int neq = Aeq.rows();
  // defulat is to have 0 linear equality constraints
  if(Aeq.size() != 0)
  {
    //Aeq = Eigen::SparseMatrix<T>(0,n);
    assert(n == Aeq.cols());
  }

  assert(A.rows() == n);
  assert(A.cols() == n);

  // number of known rows
  int kr = known.size();

  assert(kr == 0 || known.minCoeff() >= 0);
  assert(kr == 0 || known.maxCoeff() < n);
  assert(neq <= n);

  // cache known
  data.known = known;
  // get list of unknown indices
  data.unknown.resize(n-kr);
  std::vector<bool> unknown_mask;
  unknown_mask.resize(n,true);
  for(int i = 0;i<kr;i++)
  {
    unknown_mask[known(i)] = false;
  }
  int u = 0;
  for(int i = 0;i<n;i++)
  {
    if(unknown_mask[i])
    {
      data.unknown(u) = i;
      u++;
    }
  }
  // get list of lagrange multiplier indices
  data.lagrange.resize(neq);
  for(int i = 0;i<neq;i++)
  {
    data.lagrange(i) = n + i;
  }
  // cache unknown followed by lagrange indices
  data.unknown_lagrange.resize(data.unknown.size()+data.lagrange.size());
  data.unknown_lagrange << data.unknown, data.lagrange;

  Eigen::SparseMatrix<T> Auu;
  igl::slice(A,data.unknown,data.unknown,Auu);

  // determine if A(unknown,unknown) is symmetric and/or positive definite
  data.Auu_sym = igl::is_symmetric(Auu);
  // Positive definiteness is *not* determined, rather it is given as a
  // parameter
  data.Auu_pd = pd;

  // Append lagrange multiplier quadratic terms
  SparseMatrix<T> new_A;
  Eigen::Matrix<int,Eigen::Dynamic,1> AI;
  Eigen::Matrix<int,Eigen::Dynamic,1> AJ;
  Eigen::Matrix<T,Eigen::Dynamic,1> AV;
  igl::find(A,AI,AJ,AV);
  Eigen::Matrix<int,Eigen::Dynamic,1> AeqI(0,1);
  Eigen::Matrix<int,Eigen::Dynamic,1> AeqJ(0,1);
  Eigen::Matrix<T,Eigen::Dynamic,1> AeqV(0,1);
  if(neq > 0)
  {
    igl::find(Aeq,AeqI,AeqJ,AeqV);
  }
  Eigen::Matrix<int,Eigen::Dynamic,1> new_AI(AV.size()+AeqV.size()*2);
  Eigen::Matrix<int,Eigen::Dynamic,1> new_AJ(AV.size()+AeqV.size()*2);
  Eigen::Matrix<T,Eigen::Dynamic,1>   new_AV(AV.size()+AeqV.size()*2);
  new_AI << AI, (AeqI.array()+n).matrix(), AeqJ;
  new_AJ << AJ, AeqJ, (AeqI.array()+n).matrix();
  new_AV << AV, AeqV, AeqV;
  //new_AI.block(0,0,n,1) = AI;
  //new_AJ.block(0,0,n,1) = AJ;
  //new_AV.block(0,0,n,1) = AV;
  //new_AI.block(n,0,neq,1) = AeqI+n;
  //new_AJ.block(n,0,neq,1) = AeqJ;
  //new_AV.block(n,0,neq,1) = AeqV;
  //new_AI.block(n+neq,0,neq,1) = AeqJ;
  //new_AJ.block(n+neq,0,neq,1) = AeqI+n;
  //new_AV.block(n+neq,0,neq,1) = AeqV;
  igl::sparse(new_AI,new_AJ,new_AV,n+neq,n+neq,new_A);

  // precompute RHS builders
  Eigen::SparseMatrix<T> Aulk;
  igl::slice(new_A,data.unknown_lagrange,data.known,Aulk);
  Eigen::SparseMatrix<T> Akul;
  igl::slice(new_A,data.known,data.unknown_lagrange,Akul);

  //// This doesn't work!!!
  //data.preY = Aulk + Akul.transpose();
  Eigen::SparseMatrix<T> AkulT = Akul.transpose();
  //// Resize preY before assigning
  //data.preY.resize(data.unknown_lagrange.size(),data.known.size());
  data.preY = Aulk + AkulT;

  // Create factorization
  if(data.Auu_sym && data.Auu_pd)
  {
    Eigen::SparseMatrix<T> Aequ(0,0);
    if(neq>0)
    {
      Eigen::Matrix<int,Eigen::Dynamic,1> Aeq_rows;
      Aeq_rows.resize(neq);
      for(int i = 0;i<neq;i++)
      {
        Aeq_rows(i) = i;
      }
      igl::slice(Aeq,Aeq_rows,data.unknown,Aequ);
    }
    // Get transpose of Aequ
    Eigen::SparseMatrix<T> AequT = Aequ.transpose();
    // Compute LU decomposition
    bool lu_success = igl::lu_lagrange(Auu,AequT,data.L,data.U);
    if(!lu_success)
    {
      return false;
    }
  }else
  {
    Eigen::SparseMatrix<T> NA;
    igl::slice(new_A,data.unknown_lagrange,data.unknown_lagrange,NA);
    assert(false);
    return false;
  }
  return true;
}


template <typename T>
bool igl::min_quad_with_fixed_solve(
  const igl::min_quad_with_fixed_data<T> & data,
  const Eigen::Matrix<T,Eigen::Dynamic,1> & B,
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Y,
  const Eigen::Matrix<T,Eigen::Dynamic,1> & Beq,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Z)
{
  // number of known rows
  int kr = data.known.size();
  if(kr!=0)
  {
    assert(kr == Y.rows());
  }
  // number of columns to solve
  int cols = Y.cols();

  // number of lagrange multipliers aka linear equality constraints
  int neq = data.lagrange.size();

  if(neq != 0)
  {
  }

  // append lagrange multiplier rhs's
  Eigen::Matrix<T,Eigen::Dynamic,1> BBeq(B.size() + Beq.size());
  BBeq << B, (Beq*-2.0);

  // Build right hand side
  Eigen::Matrix<T,Eigen::Dynamic,1> BBequl;
  igl::slice(BBeq,data.unknown_lagrange,BBequl);
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> BBequlcols;
  igl::repmat(BBequl,1,cols,BBequlcols);
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> NB;
  NB = data.preY * Y + BBequlcols;

  // resize output
  Z.resize(data.n,cols);

  // Set known values
  for(int i = 0;i < kr;i++)
  {
    for(int j = 0;j < cols;j++)
    {
      Z(data.known(i),j) = Y(i,j);
    }
  }
  data.L.template triangularView<Lower>().solveInPlace(NB);
  data.U.template triangularView<Upper>().solveInPlace(NB);
  // Now NB contains sol/-0.5
  NB *= -0.5;
  // Now NB contains solution
  // Place solution in Z
  for(int i = 0;i<(NB.rows()-neq);i++)
  {
    for(int j = 0;j<NB.cols();j++)
    {
      Z(data.unknown_lagrange(i),j) = NB(i,j);
    }
  }
  return true;
}
#endif
