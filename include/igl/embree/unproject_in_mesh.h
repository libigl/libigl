#ifndef IGL_UNPROJECT_IN_MESH
#define IGL_UNPROJECT_IN_MESH
#include <igl/igl_inline.h>
#include <Eigen/Core>

#include <vector>
#include "Hit.h"

namespace igl
{
  // Forward define
  template <
    typename Scalar,
    typename Index>
  class EmbreeIntersector;
  // Unproject a screen location (using current opengl viewport, projection, and
  // model view) to a 3D position 
  //
  // Inputs:
  //    x  x-coordinate of mouse location
  //    y  y-coordinate of mouse location
  //    ei  EmbreeIntersector containing (V,F)
  // Outputs:
  //    obj  3d unprojected mouse point in mesh
  // Returns number of hits
  //
  template <
    typename Scalar,
    typename Index,
    typename Derivedobj>
  int unproject_in_mesh(
    const int x,
    const int y,
    const igl::EmbreeIntersector<Scalar,Index> & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj);

  template <
    typename Scalar,
    typename Index,
    typename Derivedobj>
  int unproject_in_mesh(
    const int x,
    const int y,
    const igl::EmbreeIntersector<Scalar,Index> & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj,
    std::vector<igl::Hit > & hits);
}
#ifdef IGL_HEADER_ONLY
#  include "unproject_in_mesh.cpp"
#endif
#endif
