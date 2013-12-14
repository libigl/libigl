#ifndef IGL_PARTITION_H
#define IGL_PARTITION_H
#include "igl_inline.h"
#include <Eigen/Dense>

namespace igl
{
  // PARTITION partition vertices into groups based on each
  // vertex's vector: vertices with similar coordinates (close in 
  // space) will be put in the same group.
  //
  // Inputs:
  //   W  #W by dim coordinate matrix
  //   k  desired number of groups default is dim
  // Output:
  //   G  #W list of group indices (1 to k) for each vertex, such that vertex i 
  //     is assigned to group G(i)
  //   S  k  list of seed vertices
  //   D  #W list of squared distances for each vertex to it's corresponding
  //     closest seed
  IGL_INLINE void partition(
    const Eigen::MatrixXd & W,
    const int k,
    Eigen::Matrix<int,Eigen::Dynamic,1> & G,
    Eigen::Matrix<int,Eigen::Dynamic,1> & S,
    Eigen::Matrix<double,Eigen::Dynamic,1> & D);
}

#ifdef IGL_HEADER_ONLY
#include "partition.cpp"
#endif
#endif
