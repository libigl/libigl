#ifndef IGL_BOUNDARY_CONDITIONS_H
#define IGL_BOUNDARY_CONDITIONS_H
#include "igl_inline.h"
#include <Eigen/Dense>

// Note: the elements field is currently unused but is left her for consistency
// with the matlab version (where it is also unused). 10/25/2012
namespace igl
{

  // Compute boundary conditions for automatic weights computation
  // Inputs:
  //   V  #V by dim list of domain vertices
  //   Ele  #Ele by simplex-size list of simplex indices
  //   C  #C by dim list of handle positions
  //   P  #P by 1 list of point handle indices into C
  //   BE  #BE by 2 list of bone edge indices into C
  //   CE  #CE by 2 list of cage edge indices into *P*, unused
  // Outputs:
  //   b  #b list of boundary indices (indices into V of vertices which have
  //     known, fixed values)
  //   bc #b by #weights list of known/fixed values for boundary vertices (notice
  //     the #b != #weights in general because #b will include all the
  //     intermediary samples along each bone, etc.. The ordering of the weights
  //     corresponds to [P;BE]
  // Returns true if boundary conditions make sense
  IGL_INLINE bool boundary_conditions(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele,
    const Eigen::MatrixXd & C,
    const Eigen::VectorXi & P,
    const Eigen::MatrixXi & BE,
    const Eigen::MatrixXi & CE,
    Eigen::VectorXi & b,
    Eigen::MatrixXd & bc);
}

#ifdef IGL_HEADER_ONLY
#  include "boundary_conditions.cpp"
#endif

#endif
