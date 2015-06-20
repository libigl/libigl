// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "point_mesh_squared_distance.h"
#include "mesh_to_cgal_triangle_list.h"

template <typename Kernel>
IGL_INLINE void igl::cgal::point_mesh_squared_distance(
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & sqrD,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & C)
{
  using namespace std;
  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
  Tree tree;
  vector<Triangle_3> T;
  point_mesh_squared_distance_precompute(V,F,tree,T);
  return point_mesh_squared_distance(P,tree,T,sqrD,I,C);
}

template <typename Kernel>
IGL_INLINE void igl::cgal::point_mesh_squared_distance_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > & tree,
  std::vector<CGAL::Triangle_3<Kernel> > & T)
{
  using namespace std;

  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef CGAL::Point_3<Kernel> Point_3; 
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

  // Must be 3D
  assert(V.cols() == 3);
  // Must be triangles
  assert(F.cols() == 3);

  // WTF ALERT!!!! 
  //
  // There's a bug in clang probably or at least in cgal. Without calling this
  // line (I guess invoking some compilation from <vector>), clang will vomit
  // errors inside CGAL.
  //
  // http://stackoverflow.com/questions/27748442/is-clangs-c11-support-reliable
  T.reserve(0);

  // Make list of cgal triangles
  mesh_to_cgal_triangle_list(V,F,T);
  tree.clear();
  tree.insert(T.begin(),T.end());
  tree.accelerate_distance_queries();
  // accelerate_distance_queries doesn't seem actually to do _all_ of the
  // precomputation. the tree (despite being const) will still do more
  // precomputation and reorganizing on the first call of `closest_point` or
  // `closest_point_and_primitive`. Therefor, call it once here.
  tree.closest_point_and_primitive(Point_3(0,0,0));
}

template <typename Kernel>
IGL_INLINE void igl::cgal::point_mesh_squared_distance(
  const Eigen::MatrixXd & P,
  const CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > & tree,
  const std::vector<CGAL::Triangle_3<Kernel> > & T,
  Eigen::VectorXd & sqrD,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & C)
{
  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
  typedef CGAL::Point_3<Kernel>    Point_3;
  assert(P.cols() == 3);
  const int n = P.rows();
  sqrD.resize(n,1);
  I.resize(n,1);
  C.resize(n,P.cols());
  for(int p = 0;p < n;p++)
  {
    Point_3 query(P(p,0),P(p,1),P(p,2));
    // Find closest point and primitive id
    Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
    Point_3 closest_point = pp.first;
    C(p,0) = CGAL::to_double(closest_point[0]);
    C(p,1) = CGAL::to_double(closest_point[1]);
    C(p,2) = CGAL::to_double(closest_point[2]);
    sqrD(p) = CGAL::to_double((closest_point-query).squared_length());
    I(p) = pp.second - T.begin();
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::cgal::point_mesh_squared_distance_precompute<CGAL::Epick>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Epick, CGAL::AABB_triangle_primitive<CGAL::Epick, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >::iterator, CGAL::Boolean_tag<false> > > >&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >&);
template void igl::cgal::point_mesh_squared_distance<CGAL::Epeck>( const Eigen::MatrixXd & P, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::VectorXd & sqrD, Eigen::VectorXi & I, Eigen::MatrixXd & C);
template void igl::cgal::point_mesh_squared_distance<CGAL::Epick>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
template void igl::cgal::point_mesh_squared_distance<CGAL::Simple_cartesian<double> >(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
#endif
