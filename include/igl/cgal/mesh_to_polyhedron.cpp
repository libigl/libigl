// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mesh_to_polyhedron.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
template <typename Polyhedron>
IGL_INLINE bool igl::mesh_to_polyhedron(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Polyhedron & poly)
{
  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  // Postcondition: hds is a valid polyhedral surface.
  CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(poly.hds());
  B.begin_surface(V.rows(),F.rows());
  typedef typename HalfedgeDS::Vertex   Vertex;
  typedef typename Vertex::Point Point;
  assert(V.cols() == 3 && "V must be #V by 3");
  for(int v = 0;v<V.rows();v++)
  {
    B.add_vertex(Point(V(v,0),V(v,1),V(v,2)));
  }
  assert(F.cols() == 3 && "F must be #F by 3");
  for(int f=0;f<F.rows();f++)
  {
    B.begin_facet();
    for(int c = 0;c<3;c++)
    {
      B.add_vertex_to_facet(F(f,c));
    }
    B.end_facet();
  }
  if(B.error())
  {
    B.rollback();
    return false;
  }
  B.end_surface();
  return poly.is_valid();
}
