// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "bbw.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/speye.h>
#include <igl/slice_into.h>
#include <igl/min_quad_with_fixed.h>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

#include <iostream>
#include <cstdio>


igl::bbw::BBWData::BBWData():
  partition_unity(false),
  W0(),
#ifndef IGL_NO_MOSEK
  mosek_data(),
#endif
  active_set_params(),
  qp_solver(QP_SOLVER_IGL_ACTIVE_SET),
  verbosity(0)
{
  // We know that the Bilaplacian is positive semi-definite
  active_set_params.Auu_pd = true;
}

void igl::bbw::BBWData::print()
{
  using namespace std;
  cout<<"partition_unity: "<<partition_unity<<endl;
  cout<<"W0=["<<endl<<W0<<endl<<"];"<<endl;
  cout<<"qp_solver: "<<QPSolverNames[qp_solver]<<endl;
}


template <
  typename DerivedV, 
  typename DerivedEle, 
  typename Derivedb,
  typename Derivedbc, 
  typename DerivedW>
IGL_INLINE bool igl::bbw::bbw(
  const Eigen::PlainObjectBase<DerivedV> & V, 
  const Eigen::PlainObjectBase<DerivedEle> & Ele, 
  const Eigen::PlainObjectBase<Derivedb> & b, 
  const Eigen::PlainObjectBase<Derivedbc> & bc, 
  igl::bbw::BBWData & data,
  Eigen::PlainObjectBase<DerivedW> & W
  )
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  // number of domain vertices
  int n = V.rows();
  // number of handles
  int m = bc.cols();

  SparseMatrix<typename DerivedW::Scalar> L;
  cotmatrix(V,Ele,L);
  MassMatrixType mmtype = MASSMATRIX_TYPE_VORONOI;
  if(Ele.cols() == 4)
  {
    mmtype = MASSMATRIX_TYPE_BARYCENTRIC;
  }
  SparseMatrix<typename DerivedW::Scalar> M;
  SparseMatrix<typename DerivedW::Scalar> Mi;
  massmatrix(V,Ele,mmtype,M);

  invert_diag(M,Mi);

  // Biharmonic operator
  SparseMatrix<typename DerivedW::Scalar> Q = L.transpose() * Mi * L;

  W.derived().resize(n,m);
  if(data.partition_unity)
  {
    // Not yet implemented
    assert(false);
  }else
  {
    // No linear terms
    VectorXd c = VectorXd::Zero(n);
    // No linear constraints
    SparseMatrix<typename DerivedW::Scalar> A(0,n),Aeq(0,n),Aieq(0,n);
    VectorXd uc(0,1),Beq(0,1),Bieq(0,1),lc(0,1);
    // Upper and lower box constraints (Constant bounds)
    VectorXd ux = VectorXd::Ones(n);
    VectorXd lx = VectorXd::Zero(n);
    active_set_params eff_params = data.active_set_params;
    switch(data.qp_solver)
    {
      case QP_SOLVER_IGL_ACTIVE_SET:
      {
        if(data.verbosity >= 1)
        {
          cout<<"BBW: max_iter: "<<data.active_set_params.max_iter<<endl;
          cout<<"BBW: eff_max_iter: "<<eff_params.max_iter<<endl;
        }
        if(data.verbosity >= 1)
        {
          cout<<"BBW: Computing initial weights for "<<m<<" handle"<<
            (m!=1?"s":"")<<"."<<endl;
        }
        min_quad_with_fixed_data<typename DerivedW::Scalar > mqwf;
        min_quad_with_fixed_precompute(Q,b,Aeq,true,mqwf);
        min_quad_with_fixed_solve(mqwf,c,bc,Beq,W);
        // decrement
        eff_params.max_iter--;
        bool error = false;
        // Loop over handles
#pragma omp parallel for
        for(int i = 0;i<m;i++)
        {
          // Quicker exit for openmp
          if(error)
          {
            continue;
          }
          if(data.verbosity >= 1)
          {
#pragma omp critical
            cout<<"BBW: Computing weight for handle "<<i+1<<" out of "<<m<<
              "."<<endl;
          }
          VectorXd bci = bc.col(i);
          VectorXd Wi;
          // use initial guess
          Wi = W.col(i);
          SolverStatus ret = active_set(
              Q,c,b,bci,Aeq,Beq,Aieq,Bieq,lx,ux,eff_params,Wi);
          switch(ret)
          {
            case SOLVER_STATUS_CONVERGED:
              break;
            case SOLVER_STATUS_MAX_ITER:
              cerr<<"active_set: max iter without convergence."<<endl;
              break;
            case SOLVER_STATUS_ERROR:
            default:
              cerr<<"active_set error."<<endl;
              error = true;
          }
          W.col(i) = Wi;
        }
        if(error)
        {
          return false;
        }
        break;
      }
      case QP_SOLVER_MOSEK:
      {
#ifdef IGL_NO_MOSEK
        assert(false && "Use another QPSolver. Recompile without IGL_NO_MOSEK defined.");
        cerr<<"Use another QPSolver. Recompile without IGL_NO_MOSEK defined."<<endl;
        return false;
#else
        // Loop over handles
        for(int i = 0;i<m;i++)
        {
          if(data.verbosity >= 1)
          {
            cout<<"BBW: Computing weight for handle "<<i+1<<" out of "<<m<<
              "."<<endl;
          }
          VectorXd bci = bc.col(i);
          VectorXd Wi;
          // impose boundary conditions via bounds
          slice_into(bci,b,ux);
          slice_into(bci,b,lx);
          bool r = mosek_quadprog(Q,c,0,A,lc,uc,lx,ux,data.mosek_data,Wi);
          if(!r)
          {
            return false;
          }
          W.col(i) = Wi;
        }
#endif
        break;
      }
      default:
      {
        assert(false && "Unknown qp_solver");
        return false;
      }
    }
#ifndef NDEBUG
    const double min_rowsum = W.rowwise().sum().array().abs().minCoeff();
    if(min_rowsum < 0.1)
    {
      cerr<<"bbw.cpp: Warning, minimum row sum is very low. Consider more "
        "active set iterations or enforcing partition of unity."<<endl;
    }
#endif
  }

  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template bool igl::bbw::bbw<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::bbw::BBWData&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
