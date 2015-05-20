// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LSCM_H
#define IGL_LSCM_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // Compute a Least-squares conformal map parametrization following the
  // algorithm presented in: Spectral Conformal Parameterization, Patrick
  // Mullen, Yiying Tong, Pierre Alliez and Mathieu Desbrun. Input should be a
  // manifold mesh (also no unreferenced vertices) and "boundary" `b` should
  // contain at least two vertices per connected component.
  //
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of mesh faces (must be triangles)
  //   b  #b boundary indices into V
  //   bc #b by 3 list of boundary values
  // Outputs:
  //   UV #V by 2 list of 2D mesh vertex positions in UV space
  // Returns true only on solver success.
  //
  IGL_INLINE bool lscm(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& b,
    const Eigen::MatrixXd& bc,
    Eigen::MatrixXd& V_uv);
}

#ifndef IGL_STATIC_LIBRARY
#  include "lscm.cpp"
#endif

#endif
