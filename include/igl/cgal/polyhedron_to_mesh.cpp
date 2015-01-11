// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "polyhedron_to_mesh.h"
#include <CGAL/Polyhedron_3.h>

template <typename Polyhedron>
IGL_INLINE void igl::polyhedron_to_mesh(
  const Polyhedron & poly,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F)
{
  using namespace std;
  V.resize(poly.size_of_vertices(),3);
  F.resize(poly.size_of_facets(),3);
  typedef typename Polyhedron::Vertex_const_iterator Vertex_iterator;
  std::map<Vertex_iterator,size_t> vertex_to_index;
  {
    size_t v = 0;
    for(
      typename Polyhedron::Vertex_const_iterator p = poly.vertices_begin();
      p != poly.vertices_end();
      p++)
    {
      V(v,0) = p->point().x();
      V(v,1) = p->point().y();
      V(v,2) = p->point().z();
      vertex_to_index[p] = v;
      v++;
    }
  }
  {
    size_t f = 0;
    for(
      typename Polyhedron::Facet_const_iterator facet = poly.facets_begin();
      facet != poly.facets_end();
      ++facet)
    {
      typename Polyhedron::Halfedge_around_facet_const_circulator he = 
        facet->facet_begin();
      // Facets in polyhedral surfaces are at least triangles.
      assert(CGAL::circulator_size(he) == 3 && "Facets should be triangles");
      size_t c = 0;
      do {
        //// This is stooopidly slow
        // F(f,c) = std::distance(poly.vertices_begin(), he->vertex());
        F(f,c) = vertex_to_index[he->vertex()];
        c++;
      } while ( ++he != facet->facet_begin());
      f++;
    }
  }
}
