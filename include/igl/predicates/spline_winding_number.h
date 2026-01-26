#ifndef IGL_PREDICATES_SPLINE_WINDING_NUMBER_H
#define IGL_PREDICATES_SPLINE_WINDING_NUMBER_H
#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  namespace predicates
  {
    /// @param[in] q  1 by dimensions query point
    ///
    /// \see igl::cycodebase::box_cubic
    template <
      typename DerivedP, 
      typename DerivedC, 
      typename DerivedB,
      typename Derivedleaf,
      typename DerivedQ,
      typename DerivedW>
    IGL_INLINE void spline_winding_number(
      const Eigen::MatrixBase<DerivedP>& P,
      const Eigen::MatrixBase<DerivedC>& C,
      const Eigen::MatrixBase<DerivedB>& B1,
      const Eigen::MatrixBase<DerivedB>& B2,
      const Eigen::MatrixBase<Derivedleaf>& leaf,
      const Eigen::MatrixBase<DerivedQ> & Q,
      Eigen::PlainObjectBase<DerivedW>& W);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "spline_winding_number.cpp"
#endif

#endif


