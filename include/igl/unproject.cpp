// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject.h"

#include <Eigen/Dense>
#include <Eigen/LU>

template <typename Scalar>
IGL_INLINE Eigen::Matrix<Scalar,3,1> igl::unproject(
  const    Eigen::Matrix<Scalar,3,1>&  win,
  const    Eigen::Matrix<Scalar,4,4>& model,
  const    Eigen::Matrix<Scalar,4,4>& proj,
  const    Eigen::Matrix<Scalar,4,1>&  viewport)
{
  Eigen::Matrix<Scalar,4,4> Inverse = (proj * model).inverse();

  Eigen::Matrix<Scalar,4,1> tmp;
  tmp << win, 1;
  tmp(0) = (tmp(0) - viewport(0)) / viewport(2);
  tmp(1) = (tmp(1) - viewport(1)) / viewport(3);
  tmp = tmp.array() * 2.0f - 1.0f;

  Eigen::Matrix<Scalar,4,1> obj = Inverse * tmp;
  obj /= obj(3);

  return obj.head(3);
}

#ifdef IGL_STATIC_LIBRARY
template Eigen::Matrix<float, 3, 1, 0, 3, 1> igl::unproject<float>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> igl::unproject<double>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 1, 0, 4, 1> const&);
#endif
