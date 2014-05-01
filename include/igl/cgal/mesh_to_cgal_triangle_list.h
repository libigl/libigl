#ifndef IGL_MESH_TO_CGAL_TRIANGLE_LIST_H
#define IGL_MESH_TO_CGAL_TRIANGLE_LIST_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include "CGAL_includes.hpp"
namespace igl
{
  // Convert a mesh (V,F) to a list of CGAL triangles
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   T  #F list of CGAL triangles
  template <typename Kernel>
  IGL_INLINE void mesh_to_cgal_triangle_list(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    std::vector<CGAL::Triangle_3<Kernel> > & T);
}
#ifdef IGL_HEADER_ONLY
#  include "mesh_to_cgal_triangle_list.cpp"
#endif

#endif
