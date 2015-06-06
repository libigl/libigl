// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "point_mesh_squared_distance.h"
#include "AABB.h"

template <
  typename DerivedP,
  typename DerivedV,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE void igl::point_mesh_squared_distance(
  const Eigen::PlainObjectBase<DerivedP> & P,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  using namespace std;
  AABB<DerivedV,3> tree;
  tree.init(V,Ele);
  return tree.squared_distance(V,Ele,P,sqrD,I,C);
}

#ifdef IGL_STATIC_LIBRARY
#endif
