#ifndef PEAL_OUTER_HULL_LAYERS_H
#define PEAL_OUTER_HULL_LAYERS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Computes necessary generic information for boolean operations by
  // successively "pealing" off the "outer hull" of a mesh (V,F) resulting from
  // "resolving" all (self-)intersections.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  // Outputs:
  //   odd  #F list of whether facet belongs to an odd iteration peal (otherwise
  //     an even iteration peal)
  //   flip  #F list of whether a facet's orientation was flipped when facet
  //     "pealed" into its associated outer hull layer.
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedodd,
    typename Derivedflip>
  IGL_INLINE void peal_outer_hull_layers(
    const Eigen::PlainObjectBase<DerivedV > & V,
    const Eigen::PlainObjectBase<DerivedF > & F,
    Eigen::PlainObjectBase<Derivedodd > & odd,
    Eigen::PlainObjectBase<Derivedflip > & flip);
}

#ifndef IGL_STATIC_LIBRARY
#  include "peal_outer_hull_layers.cpp"
#endif
#endif
