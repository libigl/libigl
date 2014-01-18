// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "arap.h"
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/group_sum_matrix.h>
#include <igl/covariance_scatter_matrix.h>
#include <igl/speye.h>
#include <igl/mode.h>
#include <igl/slice.h>
#include <igl/arap_rhs.h>
#include <igl/repdiag.h>
#include <igl/columnize.h>
#include "fit_rotations.h"
#include <cassert>
#include <iostream>

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedb>
IGL_INLINE bool igl::arap_precomputation(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<Derivedb> & b,
  ARAPData & data)
{
  using namespace igl;
  using namespace Eigen;
  typedef typename DerivedV::Scalar Scalar;
  // number of vertices
  const int n = V.rows();
  data.n = n;
  assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
  assert((b.size() == 0 || b.minCoeff() >=0) && "b out of bounds");
  // remember b
  data.b = b;
  assert(F.cols() == 3 && "For now only triangles");
  // dimension
  const int dim = V.cols();
  assert(dim == 3 && "Only 3d supported");
  // Defaults
  data.f_ext = Eigen::MatrixXd::Zero(n,dim);

  SparseMatrix<Scalar> L;
  cotmatrix(V,F,L);

  ARAPEnergyType eff_energy = data.energy;
  if(eff_energy == ARAP_ENERGY_TYPE_DEFAULT)
  {
    switch(F.cols())
    {
      case 3:
        if(dim == 3)
        {
          eff_energy = ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
        }else
        {
          eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
        }
        break;
      case 4:
        eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
      default:
        assert(false);
    }
  }


  // Get covariance scatter matrix, when applied collects the covariance matrices
  // used to fit rotations to during optimization
  covariance_scatter_matrix(V,F,eff_energy,data.CSM);

  // Get group sum scatter matrix, when applied sums all entries of the same
  // group according to G
  SparseMatrix<double> G_sum;
  if(data.G.size() == 0)
  {
    speye(n,G_sum);
  }else
  {
    // groups are defined per vertex, convert to per face using mode
    Eigen::Matrix<int,Eigen::Dynamic,1> GG;
    if(eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
    {
      MatrixXi GF(F.rows(),F.cols());
      for(int j = 0;j<F.cols();j++)
      {
        Matrix<int,Eigen::Dynamic,1> GFj;
        slice(data.G,F.col(j),GFj);
        GF.col(j) = GFj;
      }
      mode<int>(GF,2,GG);
    }else
    {
      GG=data.G;
    }
    //printf("group_sum_matrix()\n");
    group_sum_matrix(GG,G_sum);
  }
  SparseMatrix<double> G_sum_dim;
  repdiag(G_sum,dim,G_sum_dim);
  data.CSM = (G_sum_dim * data.CSM).eval();

  arap_rhs(V,F,eff_energy,data.K);

  SparseMatrix<double> Q = (-0.5*L).eval();
  return min_quad_with_fixed_precompute(
    Q,b,SparseMatrix<double>(),true,data.solver_data);
}

template <
  typename Derivedbc,
  typename DerivedU>
IGL_INLINE bool igl::arap_solve(
  const Eigen::PlainObjectBase<Derivedbc> & bc,
  ARAPData & data,
  Eigen::PlainObjectBase<DerivedU> & U)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  assert(data.b.size() == bc.rows());
  const int dim = bc.cols();
  const int n = data.n;
  int iter = 0;
  if(U.size() == 0)
  {
    // terrible initial guess.. should at least copy input mesh
    U = MatrixXd::Zero(data.n,dim);
  }
  MatrixXd U_prev = U;
  while(iter < data.max_iter)
  {
    U_prev = U;
    // enforce boundary conditions exactly
    for(int bi = 0;bi<bc.rows();bi++)
    {
      U.row(data.b(bi)) = bc.row(bi);
    }

    MatrixXd S = data.CSM * U.replicate(dim,1);
    MatrixXd R(dim,data.CSM.rows());
#ifdef __SSE__ // fit_rotations_SSE will convert to float if necessary
    fit_rotations_SSE(S,R);
#else
    fit_rotations(S,R);
#endif
    //for(int k = 0;k<(data.CSM.rows()/dim);k++)
    //{
    //  R.block(0,dim*k,dim,dim) = MatrixXd::Identity(dim,dim);
    //}


    // distribute group rotations to vertices in each group
    MatrixXd eff_R;
    if(data.G.size() == 0)
    {
      // copy...
      eff_R = R;
    }else
    {
      eff_R.resize(dim,dim*n);
      for(int v = 0;v<n;v++)
      {
        eff_R.block(0,dim*v,dim,dim) = 
          R.block(0,dim*data.G(v),dim,dim);
      }
    }

    VectorXd Rcol;
    columnize(eff_R,n,2,Rcol);
    VectorXd Bcol = -data.K * Rcol;
    for(int c = 0;c<dim;c++)
    {
      VectorXd Uc,Bc,bcc,Beq;
      Bc = Bcol.block(c*n,0,n,1);
      bcc = bc.col(c);
      min_quad_with_fixed_solve(
        data.solver_data,
        Bc,bcc,Beq,
        Uc);
      U.col(c) = Uc;
    }
    

    iter++;
  }
  return true;
}

#ifndef IGL_HEADER_ONLY
template bool igl::arap_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, igl::ARAPData&);
template bool igl::arap_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::ARAPData&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
