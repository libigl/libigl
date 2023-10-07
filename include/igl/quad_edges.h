#ifndef IGL_QUAD_EDGES_H
#define IGL_QUAD_EDGES_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  template <
    typename DerivedQ,
    typename DerivedE >
  IGL_INLINE void quad_edges(
    const Eigen::MatrixBase<DerivedQ> & Q,
    Eigen::PlainObjectBase<DerivedE> & E);
}

#ifndef IGL_STATIC_LIBRARY
#  include "quad_edges.cpp"
#endif

#endif

