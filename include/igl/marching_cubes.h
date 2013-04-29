#ifndef IGL_MARCHINGCUBES_H
#define IGL_MARCHINGCUBES_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // marching_cubes( values, points, x_res, y_res, z_res, vertices, faces )
  //
  // performs marching cubes reconstruction on the grid defined by values, and
  // points, and generates vertices and faces
  //
  // Input:
  //  xres, yres, zres  resolutions of the grid in x,y,z dimensions
  //  values  #number_of_grid_points x 1 array -- the scalar values of an
  //    implicit function defined on the grid points (<0 in the inside of the
  //    surface, 0 on the border, >0 outside)
  //  points  #number_of_grid_points x 3 array -- 3-D positions of the grid
  //    points, ordered in x,y,z order:
  //      points[index] = the point at (x,y,z) where :
  //      x = (index % (xres -1), 
  //      y = (index / (xres-1)) %(yres-1),
  //      z = index / (xres -1) / (yres -1) ).
  //      where x,y,z index x, y, z dimensions
  //      i.e. index = x + y*xres + z*xres*yres
  // Output:
  //   vertices  #V by 3 list of mesh vertex positions
  //   faces  #F by 3 list of mesh triangle indices
  //
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void marching_cubes(
    const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &values,
    const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 3> &points,
    const unsigned x_res,
    const unsigned y_res,
    const unsigned z_res,
    Eigen::PlainObjectBase<DerivedV> &vertices,
    Eigen::PlainObjectBase<DerivedF> &faces);
}

#ifdef IGL_HEADER_ONLY
#  include "marching_cubes.cpp"
#endif

#endif
