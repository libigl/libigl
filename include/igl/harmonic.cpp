// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "harmonic.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "invert_diag.h"
#include "min_quad_with_fixed.h"
#include <Eigen/Sparse>

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedb,
  typename Derivedbc,
  typename DerivedW>
IGL_INLINE bool igl::harmonic(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<Derivedb> & b,
  const Eigen::PlainObjectBase<Derivedbc> & bc,
  const int k,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  using namespace Eigen;
  typedef typename DerivedV::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorXS;
  SparseMatrix<Scalar> L,M,Mi;
  cotmatrix(V,F,L);
  switch(F.cols())
  {
    case 3:
      massmatrix(V,F,MASSMATRIX_TYPE_VORONOI,M);
      break;
    case 4:
    default:
      massmatrix(V,F,MASSMATRIX_TYPE_BARYCENTRIC,M);
      break;
  }
  invert_diag(M,Mi);
  SparseMatrix<Scalar> Q = -L;
  for(int p = 1;p<k;p++)
  {
    Q = (Q*Mi*-L).eval();
  }
  const VectorXS B = VectorXS::Zero(V.rows(),1);
  min_quad_with_fixed_data<Scalar> data;
  min_quad_with_fixed_precompute(Q,b,SparseMatrix<Scalar>(),true,data);
  W.resize(V.rows(),bc.cols());
  for(int w = 0;w<bc.cols();w++)
  {
    const VectorXS bcw = bc.col(w);
    VectorXS Ww;
    if(!min_quad_with_fixed_solve(data,B,bcw,VectorXS(),Ww))
    {
      return false;
    }
    W.col(w) = Ww;
  }
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&); 
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
