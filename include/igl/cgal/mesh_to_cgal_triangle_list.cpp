// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mesh_to_cgal_triangle_list.h"

#include <cassert>

template <typename Kernel>
IGL_INLINE void igl::mesh_to_cgal_triangle_list(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  std::vector<CGAL::Triangle_3<Kernel> > & T)
{
  typedef CGAL::Point_3<Kernel>    Point_3;
  typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  // Must be 3D
  assert(V.cols() == 3);
  // Must be triangles
  assert(F.cols() == 3);
  T.reserve(F.rows());
  // Loop over faces
  for(int f = 0;f<(int)F.rows();f++)
  {
    T.push_back(
      Triangle_3(
        Point_3( V(F(f,0),0), V(F(f,0),1), V(F(f,0),2)),
        Point_3( V(F(f,1),0), V(F(f,1),1), V(F(f,1),2)),
        Point_3( V(F(f,2),0), V(F(f,2),1), V(F(f,2),2))));
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::mesh_to_cgal_triangle_list<CGAL::Epeck>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > >&);
template void igl::mesh_to_cgal_triangle_list<CGAL::Epick>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > >&);
#endif
