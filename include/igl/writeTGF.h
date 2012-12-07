#ifndef IGL_WRITETGF_H
#define IGL_WRITETGF_H
#include "igl_inline.h"

#include <vector>
#include <string>
#ifndef IGL_NO_EIGEN
#include <Eigen/Dense>
#endif

namespace igl
{
  // WRITETGF
  //
  // Write a graph to a .tgf file
  //
  // Input:
  //  filename  .tgf file name
  //  V  # vertices by 3 list of vertex positions
  //  E  # edges by 2 list of edge indices
  // 
  // Assumes that graph vertices are 3 dimensional
  IGL_INLINE bool writeTGF(
    const std::string tgf_filename,
    const std::vector<std::vector<double> > & C,
    const std::vector<std::vector<int> > & E);

  #ifndef IGL_NO_EIGEN
  IGL_INLINE bool writeTGF(
    const std::string tgf_filename,
    const Eigen::MatrixXd & C,
    const Eigen::MatrixXi & E);
  #endif
}

#ifdef IGL_HEADER_ONLY
#  include "writeTGF.cpp"
#endif

#endif

