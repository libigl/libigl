// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "angles.h"
#include <cassert>

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedtheta>
void igl::angles(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<Derivedtheta>& theta)
{
  theta.resize(F.rows(),F.cols());

  auto corner = [](const Eigen::PlainObjectBase<DerivedV>& x, const Eigen::PlainObjectBase<DerivedV>& y, const Eigen::PlainObjectBase<DerivedV>& z)
  {
    Eigen::RowVector3d v1 = (x-y).normalized();
    Eigen::RowVector3d v2 = (z-y).normalized();

    return acos(v1 * v2.transpose());
  };

  for(unsigned i=0; i<F.rows(); ++i)
  {
    for(unsigned j=0; j<F.cols(); ++j)
    {
      theta(i,j) = corner(
        V.row(F(i,int(j-1+F.cols())%F.cols())),
        V.row(F(i,j)),
        V.row(F(i,(j+1+F.cols())%F.cols()))
        );
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
