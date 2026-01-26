#ifndef IGL_PREDICATES_CUBIC_WINDING_NUMBER_H
#define IGL_PREDICATES_CUBIC_WINDING_NUMBER_H
#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  namespace predicates
  {
    /// @param[in] C  4 by dimensions matrix of control points for a cubic
    /// Bézier curve
    /// @param[in] q  1 by dimensions query point
    template <typename DerivedC, typename Derivedq>
    IGL_INLINE typename DerivedC::Scalar cubic_winding_number(
      const Eigen::MatrixBase<DerivedC>& C,
      const Eigen::MatrixBase<Derivedq>& q);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "cubic_winding_number.cpp"
#endif

#endif

