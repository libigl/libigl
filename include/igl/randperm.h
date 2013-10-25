#ifndef IGL_RANDPERM_H
#define IGL_RANDPERM_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Like matlab's randperm(n) but minus 1
  //
  // Inputs:
  //   n  number of elements
  // Outputs:
  //   I  n list of rand permutation of 0:n-1
  template <typename DerivedI>
  IGL_INLINE void randperm(
    const int n,
    Eigen::PlainObjectBase<DerivedI> & I);
  template <typename DerivedI>
  IGL_INLINE Eigen::PlainObjectBase<DerivedI> randperm( const int n);
}
#ifdef IGL_HEADER_ONLY
#  include "randperm.cpp"
#endif
#endif
