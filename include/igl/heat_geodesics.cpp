// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "heat_geodesics.h"
#include "grad.h"
#include "doublearea.h"
#include "cotmatrix.h"
#include "intrinsic_delaunay_cotmatrix.h"
#include "massmatrix.h"
#include "massmatrix_intrinsic.h"
#include "boundary_facets.h"
#include "unique.h"
#include "slice.h"
#include "avg_edge_length.h"


template < typename DerivedV, typename DerivedF, typename Scalar >
IGL_INLINE bool igl::heat_geodesics_precompute(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  HeatGeodesicsData<Scalar> & data)
{
  // default t value
  const Scalar h = avg_edge_length(V,F);
  const Scalar t = h*h;
  return heat_geodesics_precompute(V,F,t,data);
}

template < typename DerivedV, typename DerivedF, typename Scalar >
IGL_INLINE bool igl::heat_geodesics_precompute(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Scalar t,
  HeatGeodesicsData<Scalar> & data)
{
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorXS;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> MatrixXS;
  igl::grad(V,F,data.Grad);
  // div
  VectorXS dblA;
  assert(F.cols() == 3 && "Only triangles are supported");
  igl::doublearea(V,F,dblA);
  data.Div = -0.25*data.Grad.transpose()*
    dblA.replicate(V.cols(),1).asDiagonal();
  Eigen::SparseMatrix<Scalar> L,M;
  Eigen::Matrix<Scalar,Eigen::Dynamic,3> l_intrinsic;
  DerivedF F_intrinsic;
  if(data.use_intrinsic_delaunay)
  {
    igl::intrinsic_delaunay_cotmatrix(V,F,L,l_intrinsic,F_intrinsic);
    igl::massmatrix_intrinsic(l_intrinsic,F_intrinsic,MASSMATRIX_TYPE_DEFAULT,M);
  }else
  {
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  }
  Eigen::SparseMatrix<Scalar> Q = M - t*L;
  Eigen::MatrixXi O;
  igl::boundary_facets(F,O);
  igl::unique(O,data.b);
  {
    Eigen::SparseMatrix<Scalar> _;
    if(!igl::min_quad_with_fixed_precompute(
      Q,Eigen::VectorXi(),_,true,data.Neumann))
    {
      return false;
    }
    // Only need if there's a boundary
    if(data.b.size()>0)
    {
      if(!igl::min_quad_with_fixed_precompute(Q,data.b,_,true,data.Dirichlet))
      {
        return false;
      }
    }
    const Eigen::SparseMatrix<double> Aeq = M.diagonal().transpose().sparseView();
    if(!igl::min_quad_with_fixed_precompute(
      (-L*0.5).eval(),Eigen::VectorXi(),Aeq,true,data.Poisson))
    {
      return false;
    }
  }
  return true;
}

template < typename Scalar, typename Derivedgamma, typename DerivedD>
IGL_INLINE void igl::heat_geodesics_solve(
  const HeatGeodesicsData<Scalar> & data,
  const Eigen::MatrixBase<Derivedgamma> & gamma,
  Eigen::PlainObjectBase<DerivedD> & D)
{
  // number of mesh vertices
  const int n = data.Grad.cols();
  // Set up delta at gamma
  DerivedD u0 = DerivedD::Zero(n,1);
  for(int g = 0;g<gamma.size();g++)
  {
    u0(gamma(g)) = 1;
  }
  // Neumann solution
  DerivedD u;
  igl::min_quad_with_fixed_solve(
    data.Neumann,u0,DerivedD(),DerivedD(),u);
  if(data.b.size()>0)
  {
    // Average Dirichelt and Neumann solutions
    DerivedD uD;
    igl::min_quad_with_fixed_solve(
      data.Dirichlet,u0,DerivedD::Zero(data.b.size()),DerivedD(),uD);
    u += uD;
    u *= 0.5;
  }
  DerivedD grad_u = data.Grad*u;
  const int m = data.Grad.rows()/3;
  for(int i = 0;i<m;i++)
  {
    const Scalar norm = sqrt(
      grad_u(0*m+i)*grad_u(0*m+i)+
      grad_u(1*m+i)*grad_u(1*m+i)+
      grad_u(2*m+i)*grad_u(2*m+i));
    if(norm == 0)
    {
      grad_u(0*m+i) = 0;
      grad_u(1*m+i) = 0;
      grad_u(2*m+i) = 0;
    }else
    {
      grad_u(0*m+i) /= norm;
      grad_u(1*m+i) /= norm;
      grad_u(2*m+i) /= norm;
    }
  }
  const DerivedD div_X = -data.Div*grad_u;
  const DerivedD Beq = (DerivedD(1,1)<<0).finished();
  igl::min_quad_with_fixed_solve(data.Poisson,-2.0*div_X,DerivedD(),Beq,D);
  DerivedD Dgamma;
  igl::slice(D,gamma,Dgamma);
  D.array() -= Dgamma.mean();
  if(D.mean() < 0)
  {
    D = -D;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::heat_geodesics_solve<double, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::HeatGeodesicsData<double> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template bool igl::heat_geodesics_precompute<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, double, igl::HeatGeodesicsData<double>&);
template bool igl::heat_geodesics_precompute<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::HeatGeodesicsData<double>&);
#endif
