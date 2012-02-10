#ifndef IGL_READDMAT_H
#define IGL_READDMAT_H
#include "igl_inline.h"
// .dmat is a simple ascii matrix file type, defined as follows. The first line
// is always:
// <#columns> <#rows>
// Then the coefficients of the matrix are given separated by whitespace with
// columns running fastest.
//
// Example:
//   The matrix m = [1 2 3; 4 5 6];
//   corresponds to a .dmat file containing:
//   3 2
//   1 4 2 5 3 6
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Read a matrix from an ascii dmat file
  //
  // Inputs:
  //   file_name  path to .dmat file
  // Outputs:
  //   W  eigen matrix containing read-in coefficients
  // Returns true on success, false on error
  //
  IGL_INLINE bool readDMAT(const std::string file_name, Eigen::MatrixXd & W);
}

#ifdef IGL_HEADER_ONLY
#  include "readDMAT.cpp"
#endif

#endif
