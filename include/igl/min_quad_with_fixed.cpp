#include "min_quad_with_fixed.h"

#include "slice.h"
#include "is_symmetric.h"
#include "find.h"
#include "sparse.h"
#include "repmat.h"
#include "lu_lagrange.h"
#include "full.h"
#include "matlab_format.h"
#include "EPS.h"
#include "cat.h"

//#include <Eigen/SparseExtra>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <cassert>
#include <cstdio>
#include <iostream>

template <typename T, typename Derivedknown>
IGL_INLINE bool igl::min_quad_with_fixed_precompute(
  const Eigen::SparseMatrix<T>& A,
  const Eigen::PlainObjectBase<Derivedknown> & known,
  const Eigen::SparseMatrix<T>& Aeq,
  const bool pd,
  min_quad_with_fixed_data<T> & data
  )
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
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

  // Positive definiteness is *not* determined, rather it is given as a
  // parameter
  data.Auu_pd = pd;
  if(data.Auu_pd)
  {
    // PD implies symmetric
    data.Auu_sym = true;
    assert(igl::is_symmetric(Auu,igl::EPS<double>()));
  }else
  {
    // determine if A(unknown,unknown) is symmetric and/or positive definite
    data.Auu_sym = igl::is_symmetric(Auu,igl::EPS<double>());
  }

  // Append lagrange multiplier quadratic terms
  SparseMatrix<T> new_A;
  SparseMatrix<T> AeqT = Aeq.transpose();
  SparseMatrix<T> Z(neq,neq);
  // This is a bit slower. But why isn't cat fast?
  new_A = 
    cat(1,
      cat(2,   A, AeqT ),
      cat(2, Aeq,    Z ));
  //cout<<matlab_format(new_A,"new_A")<<endl;
  //Matrix<int,Dynamic,1> AI;
  //Matrix<int,Dynamic,1> AJ;
  //Matrix<T,Dynamic,1> AV;
  //// Slow
  //igl::find(A,AI,AJ,AV);
  //Matrix<int,Dynamic,1> AeqI(0,1);
  //Matrix<int,Dynamic,1> AeqJ(0,1);
  //Matrix<T,Dynamic,1> AeqV(0,1);
  //if(neq > 0)
  //{
  //  // Slow
  //  igl::find(Aeq,AeqI,AeqJ,AeqV);
  //}
  //Matrix<int,Dynamic,1> new_AI(AV.size()+AeqV.size()*2);
  //Matrix<int,Dynamic,1> new_AJ(AV.size()+AeqV.size()*2);
  //Matrix<T,Dynamic,1>   new_AV(AV.size()+AeqV.size()*2);
  //new_AI << AI, (AeqI.array()+n).matrix(), AeqJ;
  //new_AJ << AJ, AeqJ, (AeqI.array()+n).matrix();
  //new_AV << AV, AeqV, AeqV;
  //// Slow
  //igl::sparse(new_AI,new_AJ,new_AV,n+neq,n+neq,new_A);
  //cout<<matlab_format(new_A,"new_A")<<endl;

  // precompute RHS builders
  if(kr > 0)
  {
    SparseMatrix<T> Aulk,Akul;
    // Slow
    igl::slice(new_A,data.unknown_lagrange,data.known,Aulk);

    //// This doesn't work!!!
    //data.preY = Aulk + Akul.transpose();
    // Slow
    if(data.Auu_sym)
    {
      data.preY = Aulk*2;
    }else
    {
      igl::slice(new_A,data.known,data.unknown_lagrange,Akul);
      SparseMatrix<T> AkulT = Akul.transpose();
      data.preY = Aulk + AkulT;
    }
  }else
  {
    data.preY.resize(data.unknown_lagrange.size(),0);
  }

  // Positive definite and no equality constraints (Postive definiteness
  // implies symmetric)
  if(data.Auu_pd && neq == 0)
  {
    data.llt.compute(Auu);
    switch(data.llt.info())
    {
      case Eigen::Success:
        break;
      case Eigen::NumericalIssue:
        cerr<<"Error: Numerical issue."<<endl;
        return false;
      default:
        cerr<<"Error: Other."<<endl;
        return false;
    }
    data.solver_type = igl::min_quad_with_fixed_data<T>::LLT; 
  }else
  {
    // Either not PD or there are equality constraints
    SparseMatrix<T> NA;
    igl::slice(new_A,data.unknown_lagrange,data.unknown_lagrange,NA);
    data.NA = NA;
    // Ideally we'd use LDLT but Eigen doesn't support positive semi-definite
    // matrices:
    // http://forum.kde.org/viewtopic.php?f=74&t=106962&p=291990#p291990
    if(data.Auu_sym && false)
    {
      data.ldlt.compute(NA);
      switch(data.ldlt.info())
      {
        case Eigen::Success:
          break;
        case Eigen::NumericalIssue:
          cerr<<"Error: Numerical issue."<<endl;
          return false;
        default:
          cerr<<"Error: Other."<<endl;
          return false;
      }
      data.solver_type = igl::min_quad_with_fixed_data<T>::LDLT; 
    }else
    {
      // Resort to LU
      // Bottleneck >1/2
      data.lu.compute(NA); 
      //std::cout<<"NA=["<<std::endl<<NA<<std::endl<<"];"<<std::endl;
      switch(data.lu.info())
      {
        case Eigen::Success:
          break;
        case Eigen::NumericalIssue:
          cerr<<"Error: Numerical issue."<<endl;
          return false;
        case Eigen::InvalidInput:
          cerr<<"Error: Invalid Input."<<endl;
          return false;
        default:
          cerr<<"Error: Other."<<endl;
          return false;
      }
      data.solver_type = igl::min_quad_with_fixed_data<T>::LU; 
    }
  }

  return true;
}


template <
  typename T,
  typename DerivedB, 
  typename DerivedY,
  typename DerivedBeq, 
  typename DerivedZ,
  typename Derivedsol>
IGL_INLINE bool igl::min_quad_with_fixed_solve(
  const min_quad_with_fixed_data<T> & data,
  const Eigen::PlainObjectBase<DerivedB> & B,
  const Eigen::PlainObjectBase<DerivedY> & Y,
  const Eigen::PlainObjectBase<DerivedBeq> & Beq,
  Eigen::PlainObjectBase<DerivedZ> & Z,
  Eigen::PlainObjectBase<Derivedsol> & sol)
{
  using namespace std;
  // number of known rows
  int kr = data.known.size();
  if(kr!=0)
  {
    assert(kr == Y.rows());
  }
  // number of columns to solve
  int cols = Y.cols();
  assert(B.cols() == 1);
  assert(Beq.size() == 0 || Beq.cols() == 1);

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

  //cout<<matlab_format(NB,"NB")<<endl;
  switch(data.solver_type)
  {
    case igl::min_quad_with_fixed_data<T>::LLT:
      sol = data.llt.solve(NB);
      break;
    case igl::min_quad_with_fixed_data<T>::LDLT:
      sol = data.ldlt.solve(NB);
      break;
    case igl::min_quad_with_fixed_data<T>::LU:
      // Not a bottleneck
      sol = data.lu.solve(NB);
      break;
    default:
      cerr<<"Error: invalid solver type"<<endl;
      return false;
  }
  //std::cout<<"sol=["<<std::endl<<sol<<std::endl<<"];"<<std::endl;
  // Now sol contains sol/-0.5
  sol *= -0.5;
  // Now sol contains solution
  // Place solution in Z
  for(int i = 0;i<(sol.rows()-neq);i++)
  {
    for(int j = 0;j<sol.cols();j++)
    {
      Z(data.unknown_lagrange(i),j) = sol(i,j);
    }
  }
  return true;
}

template <
  typename T,
  typename DerivedB, 
  typename DerivedY,
  typename DerivedBeq, 
  typename DerivedZ>
IGL_INLINE bool igl::min_quad_with_fixed_solve(
  const min_quad_with_fixed_data<T> & data,
  const Eigen::PlainObjectBase<DerivedB> & B,
  const Eigen::PlainObjectBase<DerivedY> & Y,
  const Eigen::PlainObjectBase<DerivedBeq> & Beq,
  Eigen::PlainObjectBase<DerivedZ> & Z)
{
  Eigen::PlainObjectBase<DerivedZ> sol;
  return min_quad_with_fixed_solve(data,B,Y,Beq,Z,sol);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template bool igl::min_quad_with_fixed_solve<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::min_quad_with_fixed_data<double> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template bool igl::min_quad_with_fixed_precompute<double, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, bool, igl::min_quad_with_fixed_data<double>&);
#endif

