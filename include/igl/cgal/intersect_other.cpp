// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "intersect_other.h"
#include "CGAL_includes.hpp"
#include "mesh_to_cgal_triangle_list.h"

#ifndef IGL_FIRST_HIT_EXCEPTION
#define IGL_FIRST_HIT_EXCEPTION 10
#endif

IGL_INLINE void igl::cgal::intersect_other(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U,
  const Eigen::MatrixXi & G,
  const bool first_only,
  Eigen::MatrixXi & IF)
{
  using namespace std;
  using namespace Eigen;


  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  // 3D Primitives
  typedef CGAL::Point_3<Kernel>    Point_3;
  typedef CGAL::Segment_3<Kernel>  Segment_3; 
  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef CGAL::Plane_3<Kernel>    Plane_3;
  typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron_3; 
  // 2D Primitives
  typedef CGAL::Point_2<Kernel>    Point_2;
  typedef CGAL::Segment_2<Kernel>  Segment_2; 
  typedef CGAL::Triangle_2<Kernel> Triangle_2; 
  // 2D Constrained Delaunay Triangulation types
  typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
  typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
  typedef CGAL::Exact_intersections_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> 
    CDT_2;
  typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
  // Axis-align boxes for all-pairs self-intersection detection
  typedef std::vector<Triangle_3> Triangles;
  typedef typename Triangles::iterator TrianglesIterator;
  typedef typename Triangles::const_iterator TrianglesConstIterator;
  typedef 
    CGAL::Box_intersection_d::Box_with_handle_d<double,3,TrianglesIterator> 
    Box;

  Triangles TF,TG;
  // Compute and process self intersections
  mesh_to_cgal_triangle_list(V,F,TF);
  mesh_to_cgal_triangle_list(U,G,TG);
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Box_intersection_d/Chapter_main.html#Section_63.5 
  // Create the corresponding vector of bounding boxes
  std::vector<Box> F_boxes,G_boxes;
  const auto box_up = [](Triangles & T, std::vector<Box> & boxes)
  {
    boxes.reserve(T.size());
    for ( 
      TrianglesIterator tit = T.begin(); 
      tit != T.end(); 
      ++tit)
    {
      boxes.push_back(Box(tit->bbox(), tit));
    }
  };
  box_up(TF,F_boxes);
  box_up(TG,G_boxes);
  std::list<int> lIF;
  const auto cb = [&](const Box &a, const Box &b)
  {
    using namespace std;
    // index in F and T
    int fa = a.handle()-TF.begin();
    int fb = b.handle()-TG.begin();
    const Triangle_3 & A = *a.handle();
    const Triangle_3 & B = *b.handle();
    if(CGAL::do_intersect(A,B))
    {
      // There was an intersection
      lIF.push_back(fa);
      lIF.push_back(fb);
      if(first_only)
      {
        throw IGL_FIRST_HIT_EXCEPTION;
      }
    }
  };
  try{
    CGAL::box_intersection_d(
      F_boxes.begin(), F_boxes.end(),
      G_boxes.begin(), G_boxes.end(),
      cb);
  }catch(int e)
  {
    // Rethrow if not FIRST_HIT_EXCEPTION
    if(e != IGL_FIRST_HIT_EXCEPTION)
    {
      throw e;
    }
    // Otherwise just fall through
  }

  // Convert lIF to Eigen matrix
  assert(lIF.size()%2 == 0);
  IF.resize(lIF.size()/2,2);
  {
    int i=0;
    for(
      list<int>::const_iterator ifit = lIF.begin();
      ifit!=lIF.end();
      )
    {
      IF(i,0) = (*ifit);
      ifit++; 
      IF(i,1) = (*ifit);
      ifit++;
      i++;
    }
  }

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
