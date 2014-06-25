#ifndef IGL_PARSE_RHS_H
#define IGL_PARSE_RHS_H
#include <igl/igl_inline.h>
#include <mex.h>
#include <Eigen/Dense>
namespace igl
{
  // Reads in a matrix as a double
  //
  // Inputs:
  //   prhs  points to rhs argument
  // Outputs:
  //   V  M by N matrix 
  template <typename DerivedV>
  IGL_INLINE void parse_rhs_double(
    const mxArray *prhs[], 
    Eigen::PlainObjectBase<DerivedV> & V);
  // Reads in a matrix and subtracts 1
  template <typename DerivedV>
  IGL_INLINE void parse_rhs_index(
    const mxArray *prhs[], 
    Eigen::PlainObjectBase<DerivedV> & V);
};
#ifndef IGL_STATIC_LIBRARY
#  include "parse_rhs.cpp"
#endif
#endif
