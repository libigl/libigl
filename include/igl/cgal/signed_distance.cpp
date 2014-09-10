// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "signed_distance.h"
template <typename Kernel>
IGL_INLINE typename Kernel::FT igl::signed_distance_pseudonormal(
  const CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > & tree,
  const std::vector<CGAL::Triangle_3<Kernel> > & T,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  const typename Kernel::Point_3 & q)
{
  using namespace Eigen;
  using namespace std;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename CGAL::Triangle_3<Kernel> Triangle_3;
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef typename CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef typename CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef typename CGAL::AABB_tree<AABB_triangle_traits> Tree;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

  Point_and_primitive_id pp = tree.closest_point_and_primitive(q);
  Point_3 & p = pp.first;
  const auto & qp = q-p;
  const FT sqrd = qp.squared_length();
  Vector3d v(qp.x(),qp.y(),qp.z());
  const int f = pp.second - T.begin();
  const Triangle_3 & t = *pp.second;
  // barycentric coordinates
  const auto & area = [&p,&t](const int i, const int j)->FT
  {
    return sqrt(Triangle_3(p,t.vertex(i),t.vertex(j)).squared_area());
  };
  Vector3d b(area(1,2),area(2,0),area(0,1));
  b /= b.sum();
  // Determine which normal to use
  Vector3d n;
  const double epsilon = 1e-12;
  const int type = (b.array()<=epsilon).cast<int>().sum();
  switch(type)
  {
    case 2:
      // Find vertex
      for(int c = 0;c<3;c++)
      {
        if(b(c)>epsilon)
        {
          n = VN.row(F(f,c));
          break;
        }
      }
      break;
    case 1:
      // Find edge
      for(int c = 0;c<3;c++)
      {
        if(b(c)<=epsilon)
        {
          n = EN.row(EMAP(F.rows()*c+f));
          break;
        }
      }
      break;
    default:
      assert(false && "all barycentric coords zero.");
    case 0:
      n = FN.row(f);
      break;
  }
  const double s = (v.dot(n) >= 0 ? 1. : -1.);
  return s*sqrt(sqrd);
}

template <typename Kernel>
IGL_INLINE typename Kernel::FT igl::signed_distance_winding_number(
  const CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > & tree,
  const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
  const typename Kernel::Point_3 & q)
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename CGAL::Triangle_3<Kernel> Triangle_3;
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef typename CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef typename CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef typename CGAL::AABB_tree<AABB_triangle_traits> Tree;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  Point_and_primitive_id pp = tree.closest_point_and_primitive(q);
  Point_3 & p = pp.first;
  const auto & qp = q-p;
  const Vector3d eq(q.x(),q.y(),q.z()); 
  const FT sqrd = qp.squared_length();
  const double w = hier.winding_number(eq);
  const FT s = 1.-2.*w;
  return s*sqrt(sqrd);
}

#ifdef IGL_STATIC_LIBRARY
template CGAL::Epick::FT igl::signed_distance_winding_number<CGAL::Epick>(CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Epick, CGAL::AABB_triangle_primitive<CGAL::Epick, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >::iterator, CGAL::Boolean_tag<false> > > > const&, igl::WindingNumberAABB<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, CGAL::Epick::Point_3 const&);
template CGAL::Epick::FT igl::signed_distance_pseudonormal<CGAL::Epick>(CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Epick, CGAL::AABB_triangle_primitive<CGAL::Epick, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >::iterator, CGAL::Boolean_tag<false> > > > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, CGAL::Epick::Point_3 const&);
#endif
