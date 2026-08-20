// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "ray_triangle_intersection_parameters.h"
#include <Eigen/Geometry>

template <
  typename DerivedS,
  typename DerivedD,
  typename DerivedA,
  typename Scalar>
IGL_INLINE bool igl::ray_triangle_intersection_parameters(
  const Eigen::MatrixBase<DerivedS> & source,
  const Eigen::MatrixBase<DerivedD> & dir,
  const Eigen::MatrixBase<DerivedA> & A,
  const Eigen::MatrixBase<DerivedA> & B,
  const Eigen::MatrixBase<DerivedA> & C,
  Scalar & t,
  Scalar & u,
  Scalar & v)
{
  using Vector3S = Eigen::Matrix<Scalar,3,1>;
  const Vector3S O((Scalar)source(0),(Scalar)source(1),(Scalar)source(2));
  const Vector3S d((Scalar)dir(0),(Scalar)dir(1),(Scalar)dir(2));
  const Vector3S a((Scalar)A(0),(Scalar)A(1),(Scalar)A(2));
  const Vector3S b((Scalar)B(0),(Scalar)B(1),(Scalar)B(2));
  const Vector3S c((Scalar)C(0),(Scalar)C(1),(Scalar)C(2));
  const Vector3S e1 = b - a;
  const Vector3S e2 = c - a;
  const Vector3S pv = d.cross(e2);
  const Scalar det = e1.dot(pv);
  if(det == Scalar(0)){ t = u = v = Scalar(0); return false; }
  const Scalar inv = Scalar(1) / det;
  const Vector3S tv = O - a;
  u = tv.dot(pv) * inv;
  const Vector3S qv = tv.cross(e1);
  v = d.dot(qv) * inv;
  t = e2.dot(qv) * inv;
  return true;
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace rtip_types {
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;   // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;   // RowVector3d
  using V3f = Eigen::Matrix<float,3,1,0,3,1>;
  using R3f = Eigen::Matrix<float,1,3,1,1,3>;
}}
#define IGL_RTIP(S,D,A,T) \
  template bool igl::ray_triangle_intersection_parameters<S,D,A,T>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, \
    T&, T&, T&);
IGL_RTIP(rtip_types::R3, rtip_types::R3, rtip_types::R3, double)
IGL_RTIP(rtip_types::V3, rtip_types::V3, rtip_types::V3, double)
IGL_RTIP(rtip_types::V3, rtip_types::V3, rtip_types::R3, double)
IGL_RTIP(rtip_types::R3f, rtip_types::R3f, rtip_types::R3f, float)
IGL_RTIP(rtip_types::V3f, rtip_types::V3f, rtip_types::V3f, float)
IGL_RTIP(rtip_types::V3f, rtip_types::V3f, rtip_types::R3f, float)
#undef IGL_RTIP
#endif
