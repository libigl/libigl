#include "min_quad_with_fixed.h"

#include "slice.h"
#include "is_symmetric.h"
#include "find.h"
#include "sparse.h"
#include "repmat.h"
#include "lu_lagrange.h"
#include "full.h"
#include "EPS.h"

//#include <Eigen/SparseExtra>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <cassert>
#include <cstdio>
#include <iostream>

template <typename T>
IGL_INLINE bool igl::min_quad_with_fixed_precompute(
  const Eigen::SparseMatrix<T>& A,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & known,
  const Eigen::SparseMatrix<T>& Aeq,
  const bool pd,
  igl::min_quad_with_fixed_data<T> & data
  )
{
  using namespace Eigen;
  using namespace std;
  // number of rows
  int n = A.rows();
  // cache problem size
  data.n = n;

  int neq = Aeq.rows();
  // default is to have 0 linear equality constraints
  if(Aeq.size() != 0)
  {
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

  SparseMatrix<T> Auu;
  igl::slice(A,data.unknown,data.unknown,Auu);

  // determine if A(unknown,unknown) is symmetric and/or positive definite
  data.Auu_sym = igl::is_symmetric(Auu,igl::EPS<double>());
  // Positive definiteness is *not* determined, rather it is given as a
  // parameter
  data.Auu_pd = pd;

  // Append lagrange multiplier quadratic terms
  SparseMatrix<T> new_A;
  Matrix<int,Dynamic,1> AI;
  Matrix<int,Dynamic,1> AJ;
  Matrix<T,Dynamic,1> AV;
  igl::find(A,AI,AJ,AV);
  Matrix<int,Dynamic,1> AeqI(0,1);
  Matrix<int,Dynamic,1> AeqJ(0,1);
  Matrix<T,Dynamic,1> AeqV(0,1);
  if(neq > 0)
  {
    igl::find(Aeq,AeqI,AeqJ,AeqV);
  }
  Matrix<int,Dynamic,1> new_AI(AV.size()+AeqV.size()*2);
  Matrix<int,Dynamic,1> new_AJ(AV.size()+AeqV.size()*2);
  Matrix<T,Dynamic,1>   new_AV(AV.size()+AeqV.size()*2);
  new_AI << AI, (AeqI.array()+n).matrix(), AeqJ;
  new_AJ << AJ, AeqJ, (AeqI.array()+n).matrix();
  new_AV << AV, AeqV, AeqV;
  igl::sparse(new_AI,new_AJ,new_AV,n+neq,n+neq,new_A);

  // precompute RHS builders
  if(kr > 0)
  {
    SparseMatrix<T> Aulk;
    igl::slice(new_A,data.unknown_lagrange,data.known,Aulk);
    SparseMatrix<T> Akul;
    igl::slice(new_A,data.known,data.unknown_lagrange,Akul);

    //// This doesn't work!!!
    //data.preY = Aulk + Akul.transpose();
    SparseMatrix<T> AkulT = Akul.transpose();
    data.preY = Aulk + AkulT;
  }else
  {
    data.preY.resize(data.unknown_lagrange.size(),0);
  }

  // Create factorization
  //if(data.Auu_sym && data.Auu_pd)
  //{
  //  data.sparse = true;
  //  SparseMatrix<T> Aequ(0,0);
  //  if(neq>0)
  //  {
  //    Matrix<int,Dynamic,1> Aeq_rows;
  //    Aeq_rows.resize(neq);
  //    for(int i = 0;i<neq;i++)
  //    {
  //      Aeq_rows(i) = i;
  //    }
  //    igl::slice(Aeq,Aeq_rows,data.unknown,Aequ);
  //  }
  //  // Get transpose of Aequ
  //  SparseMatrix<T> AequT = Aequ.transpose();
  //  // Compute LU decomposition
  //  // TODO: This should be using SimplicialLDLT
  //  // Q: Does SimplicialLDLT work on "Hermitian indefinite matrices" the way
  //  // that matlabs ldl does?
  //  // A: Maybe not. The eigen documentation writes, "This class provides a
  //  // LDL^T Cholesky factorizations without square root of sparse matrices
  //  // that are selfadjoint and positive definite"
  //  bool lu_success = igl::lu_lagrange(Auu,AequT,data.L,data.U);
  //  if(!lu_success)
  //  {
  //    return false;
  //  }
  //}else
  //{
  SparseMatrix<T> NA;
  igl::slice(new_A,data.unknown_lagrange,data.unknown_lagrange,NA);
  if(data.Auu_sym)
  {
    SimplicialLDLT<SparseMatrix<T> > ldlt_solver;
    ldlt_solver.compute(NA);
    if(ldlt_solver.info() != Eigen::Success)
    {
      return false;
    }
    // Implementation incomplete
    assert(false);
    cerr<<"ERROR min_quad_with_fixed_precompute implementation incomplete"<<endl;
    return false;
  }else
  {
    // Build LU decomposition of NA
    data.sparse = false;
    fprintf(stderr,
      "Warning: min_quad_with_fixed_precompute() resorting to dense LU\n");
    Matrix<T,Dynamic,Dynamic> NAfull;
    igl::full(NA,NAfull);
    data.lu = 
      FullPivLU<Matrix<T,Dynamic,Dynamic> >(NAfull);
    if(!data.lu.isInvertible())
    {
      fprintf(stderr,
        "Error: min_quad_with_fixed_precompute() LU not invertible\n");
      return false;
    }
  }
  //}

  return true;
}


template <typename T, typename DerivedY, typename DerivedZ>
IGL_INLINE bool igl::min_quad_with_fixed_solve(
  const igl::min_quad_with_fixed_data<T> & data,
  const Eigen::Matrix<T,Eigen::Dynamic,1> & B,
  const Eigen::PlainObjectBase<DerivedY> & Y,
  const Eigen::Matrix<T,Eigen::Dynamic,1> & Beq,
  Eigen::PlainObjectBase<DerivedZ> & Z)
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

  // append lagrange multiplier rhs's
  Eigen::Matrix<T,Eigen::Dynamic,1> BBeq(B.size() + Beq.size());
  BBeq << B, (Beq*-2.0);

  // Build right hand side
  Eigen::Matrix<T,Eigen::Dynamic,1> BBequl;
  igl::slice(BBeq,data.unknown_lagrange,BBequl);
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> BBequlcols;
  igl::repmat(BBequl,1,cols,BBequlcols);
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> NB;
  if(kr == 0)
  {
    NB = BBequlcols;
  }else
  {
    NB = data.preY * Y + BBequlcols;
  }

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

  //std::cout<<"NB=["<<std::endl<<NB<<std::endl<<"];"<<std::endl;

  if(data.sparse)
  {
    //std::cout<<"data.LIJV=["<<std::endl;print_ijv(data.L,1);std::cout<<std::endl<<"];"<<
    //  std::endl<<"data.L=sparse(data.LIJV(:,1),data.LIJV(:,2),data.LIJV(:,3),"<<
    //  data.L.rows()<<","<<data.L.cols()<<");"<<std::endl;
    //std::cout<<"data.UIJV=["<<std::endl;print_ijv(data.U,1);std::cout<<std::endl<<"];"<<
    //  std::endl<<"data.U=sparse(data.UIJV(:,1),data.UIJV(:,2),data.UIJV(:,3),"<<
    //  data.U.rows()<<","<<data.U.cols()<<");"<<std::endl;
    data.L.template triangularView<Eigen::Lower>().solveInPlace(NB);
    data.U.template triangularView<Eigen::Upper>().solveInPlace(NB);
  }else
  {
    NB = data.lu.solve(NB);
  }
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

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template bool igl::min_quad_with_fixed_precompute<double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::SparseMatrix<double, 0, int> const&, bool, igl::min_quad_with_fixed_data<double>&);
#endif

