// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "project.h"

template <typename Scalar>
Eigen::Matrix<Scalar,3,1> igl::project(
  const    Eigen::Matrix<Scalar,3,1>&  obj,
  const    Eigen::Matrix<Scalar,4,4>& model,
  const    Eigen::Matrix<Scalar,4,4>& proj,
  const    Eigen::Matrix<Scalar,4,1>&  viewport)
{
  Eigen::Matrix<Scalar,4,1> tmp;
  tmp << obj,1;

  tmp = model * tmp;

  tmp = proj * tmp;

  tmp = tmp.array() / tmp(3);
  tmp = tmp.array() * 0.5f + 0.5f;
  tmp(0) = tmp(0) * viewport(2) + viewport(0);
  tmp(1) = tmp(1) * viewport(3) + viewport(1);

  return tmp.head(3);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template Eigen::Matrix<double, 3, 1, 0, 3, 1> igl::project<double>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 1, 0, 4, 1> const&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> igl::project<float>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&);
#endif
