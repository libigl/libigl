// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2023 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef IGL_ISOLINES_H
#define IGL_ISOLINES_H
#include "igl_inline.h"

#include <Eigen/Core>


namespace igl
{
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedS,
    typename Derivedvals,
    typename DerivediV,
    typename DerivediE,
    typename DerivedI>
  void isolines(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedS> & S,
    const Eigen::MatrixBase<Derivedvals> & vals,
    Eigen::PlainObjectBase<DerivediV> & iV,
    Eigen::PlainObjectBase<DerivediE> & iE,
    Eigen::PlainObjectBase<DerivedI> & I);
}

#ifndef IGL_STATIC_LIBRARY
#  include "isolines.cpp"
#endif

#endif
