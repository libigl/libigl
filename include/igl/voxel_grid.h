#ifndef IGL_VOXEL_GRID_H
#define IGL_VOXEL_GRID_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
namespace igl
{
  // Construct the cell center positions of a regular voxel grid (lattice) made
  // of perfectly square voxels.
  // 
  // Inputs:
  //   box  bounding box to enclose by grid
  //   s  number of cell centers on largest side (including 2*pad_count)
  //   pad_count  number of cells beyond box
  // Outputs:
  //   GV  res(0)*res(1)*res(2) by 3 list of cell center positions
  //   res  3-long list of dimension of voxel grid
  IGL_INLINE void voxel_grid(
    const Eigen::AlignedBox3d & box, 
    const int s,
    const int pad_count,
    Eigen::MatrixXd & GV,
    Eigen::RowVector3i & side);
}
#ifndef IGL_STATIC_LIBRARY
#  include "voxel_grid.cpp"
#endif
#endif
