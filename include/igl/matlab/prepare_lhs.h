#ifndef IGL_PREPARE_LHS_H
#define IGL_PREPARE_LHS_H
#include <igl/igl_inline.h>
#include <mex.h>
#include <Eigen/Dense>
namespace igl
{
  // Writes out a matrix as a double
  //
  // Inputs:
  //   prhs  points to rhs argument
  // Outputs:
  //   V  M by N matrix 
  template <typename DerivedV>
  IGL_INLINE void prepare_lhs_double(
    const Eigen::PlainObjectBase<DerivedV> & V,
    mxArray *plhs[]);
  // Casts to logical
  template <typename DerivedV>
  IGL_INLINE void prepare_lhs_logical(
    const Eigen::PlainObjectBase<DerivedV> & V,
    mxArray *plhs[]);
  // Writes out a matrix and adds 1
  template <typename DerivedV>
  IGL_INLINE void prepare_lhs_index(
    const Eigen::PlainObjectBase<DerivedV> & V,
    mxArray *plhs[]);
};
#ifndef IGL_STATIC_LIBRARY
#  include "prepare_lhs.cpp"
#endif
#endif

