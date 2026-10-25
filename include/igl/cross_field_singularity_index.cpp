// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Derek Wingard
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "cross_field_singularity_index.h"
#include "local_basis.h"
#include "internal_angles.h"
#include "triangle_triangle_adjacency.h"
#include "is_border_vertex.h"
#include "PI.h"
#include "PlainMatrix.h"
#include <cassert>
#include <cmath>
#include <vector>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedP,
  typename DerivedI>
IGL_INLINE void igl::cross_field_singularity_index(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedP> & PD1,
  const int N,
  Eigen::PlainObjectBase<DerivedI> & index)
{
  typedef typename DerivedV::Scalar Scalar;
  typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
  assert(V.cols() == 3 && "V must be #V by 3");
  assert(F.cols() == 3 && "F must be #F by 3 (triangles)");
  assert(PD1.rows() == F.rows() && PD1.cols() == 3 && "PD1 must be #F by 3");
  assert(N >= 1 && "N must be a positive symmetry degree");
  const Scalar period = Scalar(2.0*igl::PI)/Scalar(N);

  // Per-face right-handed tangent frame (B1, B2 = B3 × B1, B3 = normal)
  PlainMatrix<DerivedV> B1,B2,B3;
  igl::local_basis(V,F,B1,B2,B3);
  // Angle of the field representative in each face frame
  const Eigen::Index m = F.rows();
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> theta(m);
  for(Eigen::Index f = 0;f<m;f++)
  {
    const RowVector3S d = PD1.row(f).template cast<Scalar>();
    assert(d.norm() > 0 && "PD1 rows must be nonzero");
    theta(f) = std::atan2(d.dot(B2.row(f)),d.dot(B1.row(f)));
  }

  // Holonomy accumulator, initialised to the angle defect
  // K_v = 2π - Σ corner angles
  const Eigen::Index n = V.rows();
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> H =
    Eigen::Matrix<Scalar,Eigen::Dynamic,1>::Constant(n,Scalar(2.0*igl::PI));
  {
    PlainMatrix<DerivedV> A;
    igl::internal_angles(V,F,A);
    for(Eigen::Index f = 0;f<m;f++)
    {
      for(int j = 0;j<3;j++)
      {
        H(F(f,j)) -= A(f,j);
      }
    }
  }

  // Each interior edge is visited once (from the face with the smaller index)
  // so that the residual enters the two endpoint stars with exactly opposite
  // signs. This makes Σ_v index(v) = Nχ hold bit-for-bit, even for a residual
  // that lands on the ±π/N wrap boundary. Edge k of face f runs from F(f,k) to
  // F(f,k+1). Walking counterclockwise around F(f,k+1) leaves f through this
  // edge, so the residual of the f→g transport is added to F(f,k+1); walking
  // counterclockwise around F(f,k) leaves g through the same edge, so the
  // residual of the g→f transport (its negative) is added to F(f,k).
  PlainMatrix<DerivedF> TT,TTi;
  igl::triangle_triangle_adjacency(F,TT,TTi);
  for(Eigen::Index f = 0;f<m;f++)
  {
    for(int k = 0;k<3;k++)
    {
      // g <= f skips boundary edges (-1), degenerate self-adjacency, and the
      // second visit of every interior edge
      const Eigen::Index g = TT(f,k);
      if(g <= f)
      {
        continue;
      }
      const Eigen::Index a = F(f,k);
      const Eigen::Index b = F(f,(k+1)%3);
      assert(F(g,TTi(f,k)) == b && F(g,(TTi(f,k)+1)%3) == a &&
        "F must be consistently oriented across every interior edge");
      // Shared edge direction, tangent to both faces
      const RowVector3S t = (V.row(b) - V.row(a)).normalized();
      const Scalar af = std::atan2(t.dot(B2.row(f)),t.dot(B1.row(f)));
      const Scalar ag = std::atan2(t.dot(B2.row(g)),t.dot(B1.row(g)));
      // Transport θ_f into g's frame by unfolding about the edge and compare
      const Scalar delta = theta(g) - (theta(f) - af + ag);
      // Principal residual in [-π/N, π/N]: the part of the rotation that is
      // not explained by the N-fold symmetry (what cross_field_mismatch
      // rounds away)
      const Scalar r = delta - period*std::round(delta/period);
      H(b) += r;
      H(a) -= r;
    }
  }

  const std::vector<bool> is_border = igl::is_border_vertex(F);
  index.setZero(n,1);
  for(Eigen::Index v = 0;v<n;v++)
  {
    if(is_border[v])
    {
      continue;
    }
    // Exact arithmetic gives a multiple of 2π/N; rounding absorbs floating
    // point error (observed ≤ 1e-13 in double on 6e5-face meshes; the
    // tolerance leaves room for float while still catching a non-integer)
    const Scalar q = H(v)/period;
    const Scalar qr = std::round(q);
    assert(std::abs(q-qr) < Scalar(1e-3) &&
      "holonomy is not a multiple of 2π/N");
    index(v) = typename DerivedI::Scalar(qr);
  }
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedP,
  typename DerivedI>
IGL_INLINE void igl::cross_field_singularity_index(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedP> & PD1,
  Eigen::PlainObjectBase<DerivedI> & index)
{
  return cross_field_singularity_index(V,F,PD1,4,index);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, const int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, const int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, const int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cross_field_singularity_index<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
