#ifndef IGL_SMOOTH_CORNER_ADJACENCY_H
#define IGL_SMOOTH_CORNER_ADJACENCY_H

#include <Eigen/Core>
#include "igl_inline.h"

namespace igl
{
  // Determine the corner-to-face adjacency relationship that can be used for
  // computing crease-aware per-corner normals.
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of triangle mesh indices into rows of V
  //   corner_threshold_radians  dihedral angle considered non-smooth (in
  //     radians)
  // Outputs:
  //   CI  #CI list of face neighbors as indices into rows of F
  //   CC  3*#F+1 list of cumulative sizes so that CC(i*3+j+1) - CC(i*3+j) is
  //     the number of faces considered smoothly incident on corner at F(i,j)
  void smooth_corner_adjacency(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const double corner_threshold_radians,
    Eigen::VectorXi & CI,
    Eigen::VectorXi & CC);
  // Determine the effective corner-to-face adjacency relationship implied by a
  // set of indexed vertex positions (FV) and normals (FV) (e.g., those read in
  // from a .obj file).
  //
  // Inputs:
  //   FV  #F by 3 list of triangle mesh indices into rows of some V
  //   FN  #F by 3 list of triangle mesh indices into rows of some N
  // Outputs:
  //   CI  #CI list of face neighbors as indices into rows of F
  //   CC  3*#F+1 list of cumulative sizes so that CC(i*3+j+1) - CC(i*3+j) is
  //     the number of faces considered smoothly incident on corner at F(i,j)
  void smooth_corner_adjacency(
    const Eigen::MatrixXi & FV,
    const Eigen::MatrixXi & FN,
    Eigen::VectorXi & CI,
    Eigen::VectorXi & CC);
}

#ifndef IGL_STATIC_LIBRARY
#  include "smooth_corner_adjacency.cpp"
#endif

#endif
