// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2023 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef IGL_ISOLINES_INTRINSIC_H
#define IGL_ISOLINES_INTRINSIC_H
#include "igl_inline.h"

#include <Eigen/Core>


namespace igl
{
  template <
    typename DerivedF,
    typename DerivedS,
    typename Derivedvals,
    typename DerivediB,
    typename DerivediFI,
    typename DerivediE,
    typename DerivedI>
  void isolines_intrinsic(
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedS> & S,
    const Eigen::MatrixBase<Derivedvals> & vals,
    Eigen::PlainObjectBase<DerivediB> & iB,
    Eigen::PlainObjectBase<DerivediFI> & iFI,
    Eigen::PlainObjectBase<DerivediE> & iE,
    Eigen::PlainObjectBase<DerivedI> & I);
  template <
    typename DerivedF,
    typename DerivedS,
    typename DeriveduE,
    typename DerivedEMAP,
    typename DeriveduEC,
    typename DeriveduEE,
    typename DerivediB,
    typename DerivediFI,
    typename DerivediE>
  void isolines_intrinsic(
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedS> & S,
    const Eigen::MatrixBase<DeriveduE> & uE,
    const Eigen::MatrixBase<DerivedEMAP> & EMAP,
    const Eigen::MatrixBase<DeriveduEC> & uEC,
    const Eigen::MatrixBase<DeriveduEE> & uEE,
    const typename DerivedS::Scalar val,
    Eigen::PlainObjectBase<DerivediB> & iB,
    Eigen::PlainObjectBase<DerivediFI> & iFI,
    Eigen::PlainObjectBase<DerivediE> & iE);
}

#ifndef IGL_STATIC_LIBRARY
#  include "isolines_intrinsic.cpp"
#endif

#endif

