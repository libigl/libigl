#ifndef IGL_MESH_WITH_SKELETON_H
#define IGL_MESH_WITH_SKELETON_H
#include "../igl_inline.h"
#include <Eigen/Dense>

namespace igl
{
// Mesh the interior of a given surface with tetrahedra which are graded (tend
// to be small near the surface and large inside) and conform to the given
// handles and samplings thereof.
//
// Inputs:
//  V  #V by 3 list of mesh vertex positions
//  F  #F by 3 list of triangle indices
//  C  #C by 3 list of vertex positions
//  P  #P list of point handle indices
//  BE #BE by 2 list of bone-edge indices
//  CE #CE by 2 list of cage-edge indices
//  samples_per_bone  #samples to add per bone
// Outputs:
//  VV  #VV by 3 list of tet-mesh vertex positions
//  TT  #TT by 4 list of tetrahedra indices
//  FF  #FF by 3 list of surface triangle indices
// Returns true only on success
bool mesh_with_skeleton(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::VectorXi & /*P*/,
  const Eigen::MatrixXi & BE,
  const Eigen::MatrixXi & CE,
  const int samples_per_bone,
  Eigen::MatrixXd & VV,
  Eigen::MatrixXi & TT,
  Eigen::MatrixXi & FF);
}


#ifdef IGL_HEADER_ONLY
#  include "mesh_with_skeleton.cpp"
#endif

#endif
