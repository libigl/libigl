// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LINE_SEARCH_H
#define IGL_LINE_SEARCH_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>

namespace igl
{
  // TODO DOCUMENTATION MISSING
  IGL_INLINE double line_search(
    Eigen::MatrixXd& x,
    const Eigen::MatrixXd& d,
    double step_size,
    std::function<double(Eigen::MatrixXd&)> energy,
    double cur_energy);

}

#ifndef IGL_STATIC_LIBRARY
#  include "line_search.cpp"
#endif

#endif
