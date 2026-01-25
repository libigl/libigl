#ifndef IGL_CUBIC_H
#define IGL_CUBIC_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  template
    <
    typename  DerivedC,
    typename  DerivedP>
  IGL_INLINE void cubic(
    const Eigen::MatrixBase<DerivedC>& C,
    const typename DerivedP::Scalar & t,
    Eigen::PlainObjectBase<DerivedP>& P);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cubic.cpp"
#endif

#endif
