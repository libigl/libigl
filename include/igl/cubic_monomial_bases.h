#ifndef IGL_CUBIC_MONOMIAL_BASES_H
#define IGL_CUBIC_MONOMIAL_BASES_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  template <
    typename DerivedC,
    typename DerivedM,
    typename DerivedD,
    typename DerivedB>
  void cubic_monomial_bases(
    const Eigen::MatrixBase<DerivedC>& C,
    Eigen::PlainObjectBase<DerivedM>& M,
    Eigen::PlainObjectBase<DerivedD>& D,
    Eigen::PlainObjectBase<DerivedB>& B);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cubic_monomial_bases.cpp"
#endif

#endif
