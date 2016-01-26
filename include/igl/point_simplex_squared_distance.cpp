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
  const Eigen::PlainObjectBase<Derivedp> & p,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedEle> & Ele,
  const typename DerivedEle::Index primitive,
  Derivedsqr_d & sqr_d,
  Eigen::PlainObjectBase<Derivedc> & c)
{
  assert(p.size() == DIM);
  assert(V.cols() == DIM);
  assert(Ele.cols() <= DIM+1);
  assert(Ele.cols() <= 3 && "Only simplices up to triangles are considered");

  typedef Derivedsqr_d Scalar;
  typedef typename Eigen::Matrix<Scalar,1,DIM> RowVectorDIMS;
  sqr_d = std::numeric_limits<Scalar>::infinity();
  const auto & set_min = 
    [&sqr_d,&c](const Scalar & sqr_d_candidate, const RowVectorDIMS & c_candidate)
  {
    if(sqr_d_candidate < sqr_d)
    {
      sqr_d = sqr_d_candidate;
      c = c_candidate;
    }
  };

  // Simplex size
  const size_t ss = Ele.cols();
  // Only one element per node
  // plane unit normal
  bool inside_triangle = false;
  Scalar d_j = std::numeric_limits<Scalar>::infinity();
  RowVectorDIMS pp;
  // Only consider triangles, and non-degenerate triangles at that
  if(ss == 3 && 
      Ele(primitive,0) != Ele(primitive,1) && 
      Ele(primitive,1) != Ele(primitive,2) && 
      Ele(primitive,2) != Ele(primitive,0))
  {
    assert(DIM == 3 && "Only implemented for 3D triangles");
    typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
    // can't be const because of annoying DIM template
    RowVector3S v10(0,0,0);
    v10.head(DIM) = (V.row(Ele(primitive,1))- V.row(Ele(primitive,0)));
    RowVector3S v20(0,0,0);
    v20.head(DIM) = (V.row(Ele(primitive,2))- V.row(Ele(primitive,0)));
    const RowVectorDIMS n = (v10.cross(v20)).head(DIM);
    Scalar n_norm = n.norm();
    if(n_norm > 0)
    {
      const RowVectorDIMS un = n/n.norm();
      // vector to plane
      const RowVectorDIMS bc = 
        1./3.*
        ( V.row(Ele(primitive,0))+
          V.row(Ele(primitive,1))+
          V.row(Ele(primitive,2)));
      const auto & v = p-bc;
      // projected point on plane
      d_j = v.dot(un);
      pp = p - d_j*un;
      // determine if pp is inside triangle
      Eigen::Matrix<Scalar,1,3> b;
      barycentric_coordinates(
            pp,
            V.row(Ele(primitive,0)),
            V.row(Ele(primitive,1)),
            V.row(Ele(primitive,2)),
            b);
      inside_triangle = fabs(fabs(b(0)) + fabs(b(1)) + fabs(b(2)) - 1.) <= 1e-10;
    }
  }
  const auto & point_point_squared_distance = [&](const RowVectorDIMS & s)
  {
    const Scalar sqr_d_s = (p-s).squaredNorm();
    set_min(sqr_d_s,s);
  };
  if(inside_triangle)
  {
    // point-triangle squared distance
    const Scalar sqr_d_j = d_j*d_j;
    //cout<<"point-triangle..."<<endl;
    set_min(sqr_d_j,pp);
  }else
  {
    if(ss >= 2)
    {
      // point-segment distance
      // number of edges
      size_t ne = ss==3?3:1;
      for(size_t x = 0;x<ne;x++)
      {
        const size_t e1 = Ele(primitive,(x+1)%ss);
        const size_t e2 = Ele(primitive,(x+2)%ss);
        const RowVectorDIMS & s = V.row(e1);
        const RowVectorDIMS & d = V.row(e2);
        // Degenerate edge
        if(e1 == e2 || (s-d).squaredNorm()==0)
        {
          // only consider once
          if(e1 < e2)
          {
            point_point_squared_distance(s);
          }
          continue;
        }
        Eigen::Matrix<Scalar,1,1> sqr_d_j_x(1,1);
        Eigen::Matrix<Scalar,1,1> t(1,1);
        project_to_line_segment(p,s,d,t,sqr_d_j_x);
        const RowVectorDIMS q = s+t(0)*(d-s);
        set_min(sqr_d_j_x(0),q);
      }
    }else
    {
      // So then Ele is just a list of points...
      assert(ss == 1);
      const RowVectorDIMS & s = V.row(Ele(primitive,0));
      point_point_squared_distance(s);
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instanciation
template void igl::point_simplex_squared_distance<2, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
template void igl::point_simplex_squared_distance<3, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
#endif
