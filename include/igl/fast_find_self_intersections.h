// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#pragma once
#ifndef FAST_FIND_SELF_INTERSECTIONS_H
#define FAST_FIND_SELF_INTERSECTIONS_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  template <
    typename DerivedV1,
    typename DerivedF1,
    typename DerivedIF,
    typename DerivedEV,
    typename DerivedEE,
    typename DerivedEI>
  IGL_INLINE bool fast_find_self_intersections(
    const Eigen::MatrixBase<DerivedV1> & V1,
    const Eigen::MatrixBase<DerivedF1> & F1,
    const bool detect_only,
    const bool first_only,
    Eigen::PlainObjectBase<DerivedIF> & IF,
    Eigen::PlainObjectBase<DerivedEV> & EV,
    Eigen::PlainObjectBase<DerivedEE> & EE,
    Eigen::PlainObjectBase<DerivedEI> & EI);
};

#ifndef IGL_STATIC_LIBRARY
#  include "fast_find_self_intersections.cpp"
#endif

#endif
