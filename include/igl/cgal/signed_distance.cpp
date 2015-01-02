// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "signed_distance.h"
#include "../per_vertex_normals.h"
#include "../per_edge_normals.h"
#include "../per_face_normals.h"
#include "../get_seconds.h"
#include "point_mesh_squared_distance.h"


template <typename Kernel>
IGL_INLINE void igl::signed_distance(
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const SignedDistanceType sign_type,
  Eigen::VectorXd & S,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & N)
{
  using namespace Eigen;
  using namespace std;
  assert(V.cols() == 3 && "V should have 3d positions");
  assert(P.cols() == 3 && "P should have 3d positions");
  assert(F.cols() == 3 && "F should have triangles");

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename CGAL::Triangle_3<Kernel> Triangle_3;
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef typename CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef typename CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef typename CGAL::AABB_tree<AABB_triangle_traits> Tree;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

  // Prepare distance computation
  Tree tree;
  vector<Triangle_3 > T;
  point_mesh_squared_distance_precompute(V,F,tree,T);

  Eigen::MatrixXd FN,VN,EN;
  Eigen::MatrixXi E;
  Eigen::VectorXi EMAP;
  WindingNumberAABB<Eigen::Vector3d> hier;
  switch(sign_type)
  {
    default:
      assert(false && "Unknown SignedDistanceType");
    case SIGNED_DISTANCE_TYPE_DEFAULT:
    case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
      hier.set_mesh(V,F);
      hier.grow();
      break;
    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
      // [Bærentzen & Aanæs 2005]
      per_face_normals(V,F,FN);
      per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
      per_edge_normals(
        V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
      N.resize(P.rows(),3);
      break;
  }

  S.resize(P.rows(),1);
  I.resize(P.rows(),1);
  C.resize(P.rows(),3);
  for(int p = 0;p<P.rows();p++)
  {
    const Point_3 q(P(p,0),P(p,1),P(p,2));
    double s,sqrd;
    Point_and_primitive_id pp;
    switch(sign_type)
    {
      default:
        assert(false && "Unknown SignedDistanceType");
      case SIGNED_DISTANCE_TYPE_DEFAULT:
      case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
        signed_distance_winding_number(tree,hier,q,s,sqrd,pp);
        break;
      case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      {
        Vector3d n(0,0,0);
        signed_distance_pseudonormal(tree,T,F,FN,VN,EN,EMAP,q,s,sqrd,pp,n);
        N.row(p) = n;
        break;
      }
    }
    I(p) = pp.second - T.begin();
    S(p) = s*sqrt(sqrd);
    C(p,0) = pp.first.x();
    C(p,1) = pp.first.y();
    C(p,2) = pp.first.z();
  }
}


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
  typename CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator> 
        > 
      >::Point_and_primitive_id pp;
  typename Kernel::FT s,sqrd;
  Eigen::Vector3d n(0,0,0);
  signed_distance_pseudonormal<Kernel>(tree,T,F,FN,VN,EN,EMAP,q,s,sqrd,pp,n);
  return s*sqrt(sqrd);
}

template <typename Kernel>
IGL_INLINE void igl::signed_distance_pseudonormal(
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
  const typename Kernel::Point_3 & q,
  typename Kernel::FT & s,
  typename Kernel::FT & sqrd,
  typename CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator> 
        > 
      >::Point_and_primitive_id & pp,
   Eigen::Vector3d & n)
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

  pp = tree.closest_point_and_primitive(q);
  Point_3 & p = pp.first;
  const auto & qp = q-p;
  sqrd = qp.squared_length();
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
  s = (v.dot(n) >= 0 ? 1. : -1.);
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
  typename CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator> 
        > 
      >::Point_and_primitive_id pp;
  typename Kernel::FT s,sqrd;
  signed_distance_winding_number<Kernel>(tree,hier,q,s,sqrd,pp);
  return s*sqrt(sqrd);
}


template <typename Kernel>
IGL_INLINE void igl::signed_distance_winding_number(
  const CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > & tree,
  const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
  const typename Kernel::Point_3 & q,
  typename Kernel::FT & s,
  typename Kernel::FT & sqrd,
  typename CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator> 
        > 
      >::Point_and_primitive_id & pp)
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
  pp = tree.closest_point_and_primitive(q);
  Point_3 & p = pp.first;
  const auto & qp = q-p;
  const Vector3d eq(q.x(),q.y(),q.z()); 
  sqrd = qp.squared_length();
  const double w = hier.winding_number(eq);
  s = 1.-2.*w;
}

#ifdef IGL_STATIC_LIBRARY
// This template is necessary for the others to compile with clang
// http://stackoverflow.com/questions/27748442/is-clangs-c11-support-reliable
template void igl::signed_distance<CGAL::Epick>( const Eigen::MatrixXd & , const Eigen::MatrixXd & , const Eigen::MatrixXi & , const SignedDistanceType , Eigen::VectorXd & , Eigen::VectorXi &, Eigen::MatrixXd & , Eigen::MatrixXd & );
template CGAL::Epick::FT igl::signed_distance_winding_number<CGAL::Epick>(CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Epick, CGAL::AABB_triangle_primitive<CGAL::Epick, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >::iterator, CGAL::Boolean_tag<false> > > > const&, igl::WindingNumberAABB<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, CGAL::Epick::Point_3 const&);
template CGAL::Epick::FT igl::signed_distance_pseudonormal<CGAL::Epick>(CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Epick, CGAL::AABB_triangle_primitive<CGAL::Epick, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >::iterator, CGAL::Boolean_tag<false> > > > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, CGAL::Epick::Point_3 const&);
#endif
