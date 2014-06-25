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

IGL_INLINE bool igl::harmonic(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  const Eigen::MatrixXd & bc,
  const int k,
  Eigen::MatrixXd & W)
{
  using namespace igl;
  using namespace Eigen;
  SparseMatrix<double> L,M,Mi;
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
  SparseMatrix<double> Q = -L;
  for(int p = 1;p<k;p++)
  {
    Q = (Q*Mi*-L).eval();
  }
  const VectorXd B = VectorXd::Zero(V.rows(),1);
  min_quad_with_fixed_data<double> data;
  min_quad_with_fixed_precompute(Q,b,SparseMatrix<double>(),true,data);
  W.resize(V.rows(),bc.cols());
  for(int w = 0;w<bc.cols();w++)
  {
    const VectorXd bcw = bc.col(w);
    VectorXd Ww;
    if(!min_quad_with_fixed_solve(data,B,bcw,VectorXd(),Ww))
    {
      return false;
    }
    W.col(w) = Ww;
  }
  return true;
}
