#ifndef IGL_COMPONENTS_H
#define IGL_COMPONENTS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // Compute connected components of a graph represented by an adjacency matrix
  //
  // Inputs:
  //   A  n by n adjacency matrix
  // Outputs:
  //   C  n list of component ids
  //
  template <typename AScalar, typename DerivedC>
  IGL_INLINE void components(
    const Eigen::SparseMatrix<AScalar> & A,
    Eigen::PlainObjectBase<DerivedC> & C);
  // Ditto but for mesh faces as input
  //
  // Inputs:
  //   F  n by 3 list of triangle indices
  // Outputs:
  //   C  max(F) list of component ids
  template <typename DerivedF, typename DerivedC>
  IGL_INLINE void components(
    const Eigen::PlainObjectBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedC> & C);

}

#ifdef IGL_HEADER_ONLY
#  include "components.cpp"
#endif

#endif
