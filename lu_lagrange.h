#ifndef IGL_LU_LAGRANGE_H
#define IGL_LU_LAGRANGE_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // LU_LAGRANGE Compute a LU decomposition for a special type of
  // matrix Q that is symmetric but not positive-definite:
  // Q = [A'*A C
  //      C'   0];
  // where A'*A, or ATA, is given as a symmetric positive definite matrix and C
  // has full column-rank(?)
  //
  // [J] = lu_lagrange(ATA,C)
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   ATA   n by n square, symmetric, positive-definite system matrix, usually
  //     the quadratic coefficients corresponding to the original unknowns in a
  //     system
  //   C  n by m rectangular matrix corresponding the quadratic coefficients of
  //     the original unknowns times the lagrange multipliers enforcing linear
  //     equality constraints
  // Outputs:
  //   L  lower triangular matrix such that Q = L*U
  //   U  upper triangular matrix such that Q = L*U
  // Returns true on success, false on error
  //
  template <typename T>
  bool lu_lagrange(
    const SparseMatrix<T> & ATA,
    const SparseMatrix<T> & C,
    SparseMatrix<T> & L,
    SparseMatrix<T> & U);

}

// Implementation
// Cholesky LLT decomposition for symmetric positive definite
#include <Eigen/SparseExtra>
#include <cassert>
#include "find.h"
#include "sparse.h"

template <typename T>
bool igl::lu_lagrange(
  const SparseMatrix<T> & ATA,
  const SparseMatrix<T> & C,
  SparseMatrix<T> & L,
  SparseMatrix<T> & U)
{
  // number of unknowns
  int n = ATA.rows();
  // number of lagrange multipliers
  int m = C.cols();

  assert(ATA.cols() == n);
  if(m != 0)
  {
    assert(C.rows() == n);
  }


  // Cholesky factorization of ATA
  //// Eigen fails if you give a full view of the matrix like this:
  //Eigen::SparseLLT<SparseMatrix<T> > ATA_LLT(ATA);
  SparseMatrix<T> ATA_LT = ATA.template triangularView<Eigen::Lower>();
  Eigen::SparseLLT<SparseMatrix<T> > ATA_LLT(ATA_LT);

  Eigen::SparseMatrix<T> J = ATA_LLT.matrixL();

  //if(!ATA_LLT.succeeded())
  if(!((J*0).eval().nonZeros() == 0))
  {
    fprintf(stderr,"Error: lu_lagrange() failed to factor ATA\n");
    return false;
  }

  if(m == 0)
  {
    L = J;
    U = J.transpose();
  }else
  {
    // Construct helper matrix M
    Eigen::SparseMatrix<T> M = C;
    J.template triangularView<Eigen::Lower>().solveInPlace(M);


    // Compute cholesky factorizaiton of M'*M
    Eigen::SparseMatrix<T> MTM = M.transpose() * M;


    Eigen::SparseLLT<SparseMatrix<T> > MTM_LLT(MTM.template triangularView<Eigen::Lower>());

    Eigen::SparseMatrix<T> K = MTM_LLT.matrixL();


    //if(!MTM_LLT.succeeded())
    if(!((K*0).eval().nonZeros() == 0))
    {
      fprintf(stderr,"Error: lu_lagrange() failed to factor MTM\n");
      return false;
    }

    // assemble LU decomposition of Q
    Eigen::Matrix<int,Eigen::Dynamic,1> MI;
    Eigen::Matrix<int,Eigen::Dynamic,1> MJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> MV;
    igl::find(M,MI,MJ,MV);

    Eigen::Matrix<int,Eigen::Dynamic,1> KI;
    Eigen::Matrix<int,Eigen::Dynamic,1> KJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> KV;
    igl::find(K,KI,KJ,KV);

    Eigen::Matrix<int,Eigen::Dynamic,1> JI;
    Eigen::Matrix<int,Eigen::Dynamic,1> JJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> JV;
    igl::find(J,JI,JJ,JV);

    int nnz = JV.size()  + MV.size() + KV.size();

    Eigen::Matrix<int,Eigen::Dynamic,1> UI(nnz);
    Eigen::Matrix<int,Eigen::Dynamic,1> UJ(nnz);
    Eigen::Matrix<T,Eigen::Dynamic,1> UV(nnz);
    UI << JJ,                        MI, (KJ.array() + n).matrix();
    UJ << JI, (MJ.array() + n).matrix(), (KI.array() + n).matrix(); 
    UV << JV,                        MV,                     KV*-1;
    igl::sparse(UI,UJ,UV,U);

    Eigen::Matrix<int,Eigen::Dynamic,1> LI(nnz);
    Eigen::Matrix<int,Eigen::Dynamic,1> LJ(nnz);
    Eigen::Matrix<T,Eigen::Dynamic,1> LV(nnz);
    LI << JI, (MJ.array() + n).matrix(), (KI.array() + n).matrix();
    LJ << JJ,                        MI, (KJ.array() + n).matrix(); 
    LV << JV,                        MV,                        KV;
    igl::sparse(LI,LJ,LV,L);

    //// Only keep lower and upper parts
    //L = L.template triangularView<Lower>();
    //U = U.template triangularView<Upper>();
  }

  return true;
}

#endif
