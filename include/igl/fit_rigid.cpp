// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "fit_rigid.h"
#include "procrustes.h"

IGL_INLINE void igl::fit_rigid(
  const Eigen::MatrixXd & A,
  const Eigen::MatrixXd & B,
  Eigen::Rotation2Dd & R,
  Eigen::RowVector2d & t)
{
  using namespace Eigen;
  Matrix2d Rmat;
  procrustes(A,B,Rmat,t);
  R.fromRotationMatrix(Rmat);
}
