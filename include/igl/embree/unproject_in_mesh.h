#ifndef IGL_UNPROJECT_IN_MESH
#define IGL_UNPROJECT_IN_MESH

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace igl
{
  // Forward define
  template <
    typename PointMatrixType,
    typename FaceMatrixType,
    typename RowVector3>
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
  // Returns true only if ray through (x,y) hits (V,F) at least
  // once
  //
  template <
    typename PointMatrixType,
    typename FaceMatrixType,
    typename RowVector3,
    typename Derivedobj>
  bool unproject_in_mesh(
    const int x,
    const int y,
    const igl::EmbreeIntersector<PointMatrixType,FaceMatrixType,RowVector3> & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj);
}
#ifdef IGL_HEADER_ONLY
#  include "unproject_in_mesh.cpp"
#endif
#endif
