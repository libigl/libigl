// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Qingnan Zhou <qnzhou@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "orient2d.h"
#include "exactinit.h"
#include "../parallel_for.h"
#include <predicates.h>

namespace igl {
namespace predicates {

using REAL = IGL_PREDICATES_REAL;
#include "IGL_PREDICATES_ASSERT_SCALAR.h"

template<typename Vector2D>
IGL_INLINE Orientation orient2d(
    const Eigen::MatrixBase<Vector2D>& pa,
    const Eigen::MatrixBase<Vector2D>& pb,
    const Eigen::MatrixBase<Vector2D>& pc)
{
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector2D, 2);
  IGL_PREDICATES_ASSERT_SCALAR(Vector2D);

  using Point = Eigen::Matrix<REAL, 2, 1>;
  Point a{pa[0], pa[1]};
  Point b{pb[0], pb[1]};
  Point c{pc[0], pc[1]};

  const auto r = ::orient2d(a.data(), b.data(), c.data());

  if (r > 0) return Orientation::POSITIVE;
  else if (r < 0) return Orientation::NEGATIVE;
  else return Orientation::COLLINEAR;
}

template 
  <typename DerivedA,
   typename DerivedB,
   typename DerivedC,
   typename DerivedR>
IGL_INLINE void orient2d(
    const Eigen::MatrixBase<DerivedA>& A,
    const Eigen::MatrixBase<DerivedB>& B,
    const Eigen::MatrixBase<DerivedC>& C,
    Eigen::PlainObjectBase<DerivedR>& R)
{
  igl::predicates::exactinit();
  typedef typename DerivedR::Scalar RScalar;
  typedef typename DerivedA::Scalar Scalar;
  typedef Eigen::Matrix<Scalar, 1, 2> RowVector2S;

  // max(A.rows(),B.rows(),C.rows()) is the number of points
  const int np = std::max(
    std::max(A.rows(), B.rows()),C.rows());
  R.resize(np, 1);
  igl::parallel_for(np, [&](const int p)
    {
      // Not sure if these copies are needed
      const RowVector2S a = A.row(p % A.rows());
      const RowVector2S b = B.row(p % B.rows());
      const RowVector2S c = C.row(p % C.rows());
      // Compute the orientation
      R(p) = static_cast<RScalar>(igl::predicates::orient2d(a, b, c));
    },1000);
}

}
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_ORIENT2D(Vector) template igl::predicates::Orientation igl::predicates::orient2d<Vector>(const Eigen::MatrixBase<Vector>&, const Eigen::MatrixBase<Vector>&, const Eigen::MatrixBase<Vector>&)
#define IGL_MATRIX(T, R, C) Eigen::Matrix<T, R, C>
IGL_ORIENT2D(IGL_MATRIX(float, 1, 2));
IGL_ORIENT2D(IGL_MATRIX(float, 2, 1));
#ifndef LIBIGL_PREDICATES_USE_FLOAT
IGL_ORIENT2D(IGL_MATRIX(double, 1, 2));
IGL_ORIENT2D(IGL_MATRIX(double, 2, 1));
#endif
#undef IGL_MATRIX
#undef IGL_ORIENT2D
#endif
