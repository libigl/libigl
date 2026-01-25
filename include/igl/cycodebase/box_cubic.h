#ifndef IGL_CYCODEBASE_BOX_CUBIC_H
#define IGL_CYCODEBASE_BOX_CUBIC_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
  namespace cycodebase {
    template < 
      typename DerivedC, 
      typename DerivedB>
    IGL_INLINE void box_cubic(
      const Eigen::MatrixBase<DerivedC>& C,
      Eigen::PlainObjectBase<DerivedB>& B1,
      Eigen::PlainObjectBase<DerivedB>& B2);
    template < 
      typename DerivedP, 
      typename DerivedC, 
      typename DerivedB>
    IGL_INLINE void box_cubic(
      const Eigen::MatrixBase<DerivedP>& P,
      const Eigen::MatrixBase<DerivedC>& C,
      Eigen::PlainObjectBase<DerivedB>& B1,
      Eigen::PlainObjectBase<DerivedB>& B2);
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "box_cubic.cpp"
#endif
#endif
