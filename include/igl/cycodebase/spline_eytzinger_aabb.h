#ifndef IGL_CYCODEBASE_SPLINE_EYTZINGER_AABB_H
#define IGL_CYCODEBASE_SPLINE_EYTZINGER_AABB_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
  namespace cycodebase {
    template <
      typename DerivedP,
      typename DerivedC,
      typename DerivedB,
      typename Derivedleaf>
    void spline_eytzinger_aabb(
      const Eigen::MatrixBase<DerivedP>& P,
      const Eigen::MatrixBase<DerivedC>& C,
      Eigen::PlainObjectBase<DerivedB>& B1,
      Eigen::PlainObjectBase<DerivedB>& B2,
      Eigen::PlainObjectBase<Derivedleaf>& leaf);
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "spline_eytzinger_aabb.cpp"
#endif

#endif




