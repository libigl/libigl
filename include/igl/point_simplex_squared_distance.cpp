// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "point_simplex_squared_distance.h"
#include "project_to_line_segment.h"
#include "barycentric_coordinates.h"
#include <Eigen/Geometry>
#include <limits>
#include <cassert>



template <
  int DIM,
  typename Derivedp,
  typename DerivedV,
  typename DerivedEle,
  typename Derivedsqr_d,
  typename Derivedc>
IGL_INLINE void igl::point_simplex_squared_distance(
  const Eigen::MatrixBase<Derivedp> & p,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const typename DerivedEle::Index primitive,
  Derivedsqr_d & sqr_d,
  Eigen::PlainObjectBase<Derivedc> & c)
{
  typedef typename Derivedp::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar,1,DIM> Vector;
  typedef Vector Point;

  const auto & Dot = [](const Point & a, const Point & b)->Scalar
  {
    return a.dot(b);
  };
  // Real-time collision detection, Ericson, Chapter 5
  const auto & ClosestPtPointTriangle = 
    [&Dot](Point p, Point a, Point b, Point c)->Point 
  {
    // Check if P in vertex region outside A
    Vector ab = b - a;
    Vector ac = c - a;
    Vector ap = p - a;
    Scalar d1 = Dot(ab, ap);
    Scalar d2 = Dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a; // barycentric coordinates (1,0,0)
    // Check if P in vertex region outside B
    Vector bp = p - b;
    Scalar d3 = Dot(ab, bp);
    Scalar d4 = Dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b; // barycentric coordinates (0,1,0)
    // Check if P in edge region of AB, if so return projection of P onto AB
    Scalar vc = d1*d4 - d3*d2;
    if( a != b)
    {
      if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        Scalar v = d1 / (d1 - d3);
        return a + v * ab; // barycentric coordinates (1-v,v,0)
      }
    }
    // Check if P in vertex region outside C
    Vector cp = p - c;
    Scalar d5 = Dot(ab, cp);
    Scalar d6 = Dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c; // barycentric coordinates (0,0,1)
    // Check if P in edge region of AC, if so return projection of P onto AC
    Scalar vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
      Scalar w = d2 / (d2 - d6);
      return a + w * ac; // barycentric coordinates (1-w,0,w)
    }
    // Check if P in edge region of BC, if so return projection of P onto BC
    Scalar va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
      Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      return b + w * (c - b); // barycentric coordinates (0,1-w,w)
    }
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    Scalar denom = 1.0 / (va + vb + vc);
    Scalar v = vb * denom;
    Scalar w = vc * denom;
    return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
  };

  assert(p.size() == DIM);
  assert(V.cols() == DIM);
  assert(Ele.cols() <= DIM+1);
  assert(Ele.cols() <= 3 && "Only simplices up to triangles are considered");

  c = ClosestPtPointTriangle(
    p,
    V.row(Ele(primitive,0)),
    V.row(Ele(primitive,1%Ele.cols())),
    V.row(Ele(primitive,2%Ele.cols())));
  sqr_d = (p-c).squaredNorm();
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instanciation
template void igl::point_simplex_squared_distance<3, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::point_simplex_squared_distance<2, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
#endif
