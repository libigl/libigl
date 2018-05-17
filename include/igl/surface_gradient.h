// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SURFACE_GRADIENT_H
#define IGL_SURFACE_GRADIENT_H

#include "igl_inline.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  IGL_INLINE void surface_gradient(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                                  const Eigen::MatrixXd &F1, const Eigen::MatrixXd &F2,
                                                  Eigen::SparseMatrix<double> &D1, Eigen::SparseMatrix<double> &D2);
}
#ifndef IGL_STATIC_LIBRARY
#  include "surface_gradient.cpp"
#endif

#endif