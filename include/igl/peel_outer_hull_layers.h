#ifndef PEEL_OUTER_HULL_LAYERS_H
#define PEEL_OUTER_HULL_LAYERS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Computes necessary generic information for boolean operations by
  // successively "peeling" off the "outer hull" of a mesh (V,F) resulting from
  // "resolving" all (self-)intersections.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  //   N  #F by 3 list of per-face normals
  // Outputs:
  //   odd  #F list of whether facet belongs to an odd iteration peel (otherwise
  //     an even iteration peel)
  //   flip  #F list of whether a facet's orientation was flipped when facet
  //     "peeled" into its associated outer hull layer.
  // Returns number of peels
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename Derivedodd,
    typename Derivedflip>
  IGL_INLINE size_t peel_outer_hull_layers(
    const Eigen::PlainObjectBase<DerivedV > & V,
    const Eigen::PlainObjectBase<DerivedF > & F,
    const Eigen::PlainObjectBase<DerivedN > & N,
    Eigen::PlainObjectBase<Derivedodd > & odd,
    Eigen::PlainObjectBase<Derivedflip > & flip);
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedodd,
    typename Derivedflip>
  IGL_INLINE size_t peel_outer_hull_layers(
    const Eigen::PlainObjectBase<DerivedV > & V,
    const Eigen::PlainObjectBase<DerivedF > & F,
    Eigen::PlainObjectBase<Derivedodd > & odd,
    Eigen::PlainObjectBase<Derivedflip > & flip);
}

#ifndef IGL_STATIC_LIBRARY
#  include "peel_outer_hull_layers.cpp"
#endif
#endif
