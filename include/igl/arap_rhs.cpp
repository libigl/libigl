// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "arap_rhs.h"
#include "arap_linear_block.h"
#include "verbose.h"
#include "repdiag.h"
#include "cat.h"

IGL_INLINE void igl::arap_rhs(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const igl::ARAPEnergyType energy,
  Eigen::SparseMatrix<double>& K)
{
  using namespace igl;
  using namespace Eigen;
  // Number of dimensions
  int dim = V.cols();
  //// Number of mesh vertices
  //int n = V.rows();
  //// Number of mesh elements
  //int m = F.rows();
  //// number of rotations
  //int nr;
  switch(energy)
  {
    case ARAP_ENERGY_TYPE_SPOKES:
      //nr = n;
      break;
    case ARAP_ENERGY_TYPE_SPOKES_AND_RIMS:
      //nr = n;
      break;
    case ARAP_ENERGY_TYPE_ELEMENTS:
      //nr = m;
      break;
    default:
      fprintf(
        stderr,
        "covariance_scatter_matrix.h: Error: Unsupported arap energy %d\n",
        energy);
      return;
  }

  SparseMatrix<double> KX,KY,KZ;
  arap_linear_block(V,F,0,energy,KX);
  arap_linear_block(V,F,1,energy,KY);
  if(dim == 2)
  {
    K = cat(2,repdiag(KX,dim),repdiag(KY,dim));
  }else if(dim == 3)
  {
    arap_linear_block(V,F,2,energy,KZ);
    K = cat(2,cat(2,repdiag(KX,dim),repdiag(KY,dim)),repdiag(KZ,dim));
  }else
  {
    fprintf(
     stderr,
     "covariance_scatter_matrix.h: Error: Unsupported dimension %d\n",
     dim);
    return;
  }
  
}

