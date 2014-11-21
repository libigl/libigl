#ifndef FACET_COMPONENTS_H
#define FACET_COMPONENTS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{
  // Compute connected components of facets based on edge-edge adjacency.
  //
  // Inputs:
  //   F  #F by 3 list of triangle indices
  // Ouputs:
  //   C  #F list of connected component ids
  template <typename DerivedF, typename DerivedC>
  IGL_INLINE void facet_components(
    const Eigen::PlainObjectBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedC> & C);
  // Inputs:
  //   TT  #TT by 3 list of list of adjacency triangles (see
  //   triangle_triangle_adjacency.h)
  // Ouputs:
  //   C  #F list of connected component ids
  template <
    typename TTIndex, 
    typename DerivedC,
    typename Derivedcounts>
  IGL_INLINE void facet_components(
    const std::vector<std::vector<std::vector<TTIndex > > > & TT,
    Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<Derivedcounts> & counts);
}
#ifndef IGL_STATIC_LIBRARY
#  include "facet_components.cpp"
#endif
#endif
