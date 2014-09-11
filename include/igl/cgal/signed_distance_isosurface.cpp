// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "signed_distance_isosurface.h"
#include "point_mesh_squared_distance.h"
#include "complex_to_mesh.h"
#include "signed_distance.h"

#include <igl/per_face_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/centroid.h>
#include <igl/WindingNumberAABB.h>
#include <igl/matlab_format.h>
#include <igl/remove_unreferenced.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
// Axis-aligned bounding box tree for tet tri intersection
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <vector>

IGL_INLINE bool igl::signed_distance_isosurface(
  const Eigen::MatrixXd & IV,
  const Eigen::MatrixXi & IF,
  const double level,
  const double angle_bound,
  const double radius_bound,
  const double distance_bound,
  const SignedDistanceType sign_type,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F)
{
  using namespace std;

  // default triangulation for Surface_mesher
  typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
  // c2t3
  typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
  typedef Tr::Geom_traits GT;//Kernel
  typedef GT::Sphere_3 Sphere_3;
  typedef GT::Point_3 Point_3;
  typedef GT::FT FT;
  typedef std::function<FT (Point_3)> Function;
  typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
  typedef CGAL::Polyhedron_3<GT> Polyhedron;
  typedef GT::Kernel Kernel;
  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef typename std::vector<Triangle_3>::iterator Iterator;
  typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

  Tree tree;
  vector<Triangle_3 > T;
  point_mesh_squared_distance_precompute(IV,IF,tree,T);

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
      hier.set_mesh(IV,IF);
      hier.grow();
      break;
    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
      // [Bærentzen & Aanæs 2005]
      per_face_normals(IV,IF,FN);
      per_vertex_normals(IV,IF,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
      per_edge_normals(
        IV,IF,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
      break;
  }

  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  // defining the surface
  const auto & IVmax = IV.colwise().maxCoeff();
  const auto & IVmin = IV.colwise().minCoeff();
  const double bbd = (IVmax-IVmin).norm();
  const double r = bbd/2.;
  const auto & IVmid = 0.5*(IVmax + IVmin);
  // Supposedly the implict needs to evaluate to <0 at cmid...
  // http://doc.cgal.org/latest/Surface_mesher/classCGAL_1_1Implicit__surface__3.html
  Point_3 cmid(IVmid(0),IVmid(1),IVmid(2));
  Function fun;
  switch(sign_type)
  {
    default:
      assert(false && "Unknown SignedDistanceType");
    case SIGNED_DISTANCE_TYPE_DEFAULT:
    case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
      fun = 
        [&tree,&hier,&level](const Point_3 & q) -> FT
        {
          return signed_distance_winding_number(tree,hier,q)-level;
        };
      break;
    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      fun = [&tree,&T,&IF,&FN,&VN,&EN,&EMAP,&level](const Point_3 & q) -> FT
        {
          return 
            igl::signed_distance_pseudonormal(tree,T,IF,FN,VN,EN,EMAP,q) - 
            level;
        };
      break;
  }
    //[&tree,&hier,&T,&IF,&FN,&VN,&EN,&EMAP,&level](const Point_3 & q) ->FT
    //{
    //  const FT p = signed_distance_pseudonormal(tree,T,IF,FN,VN,EN,EMAP,q);
    //  const FT w = signed_distance_winding_number(tree,hier,q);
    //  if(w*p < 0 && (fabs(w) > 0.1 || fabs(p) > 0.1))
    //  {
    //    cout<<"q=["<<q.x()<<","<<q.y()<<","<<q.z()<<"];"<<endl;
    //    cout<<matlab_format(n.transpose().eval(),"n")<<endl;
    //    cout<<matlab_format(b.transpose().eval(),"b")<<endl;
    //    cout<<"Sig difference: "<<type<<endl;
    //    cout<<"w: "<<w<<endl;
    //    cout<<"p: "<<p<<endl;
    //    exit(1);
    //  }
    //  return w;
    //},
  Sphere_3 bounding_sphere(cmid, (r+level)*(r+level));
  Surface_3 surface(fun,bounding_sphere);
  CGAL::Surface_mesh_default_criteria_3<Tr> 
    criteria(angle_bound,radius_bound,distance_bound);
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  // complex to (V,F)
  return igl::complex_to_mesh(c2t3,V,F);
}
