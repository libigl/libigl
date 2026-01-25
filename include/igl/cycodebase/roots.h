#ifndef IGL_CYCODEBASE_ROOTS_H
#define IGL_CYCODEBASE_ROOTS_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
  namespace cycodebase {
    template <
      typename Derivedcoef,
      typename DerivedR>
    IGL_INLINE int roots(
        const Eigen::MatrixBase<Derivedcoef>& coef,
        const typename Derivedcoef::Scalar xmin,
        const typename Derivedcoef::Scalar xmax,
        Eigen::PlainObjectBase<DerivedR>& R);
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "roots.cpp"
#endif

#endif
