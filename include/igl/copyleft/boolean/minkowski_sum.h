#ifndef IGL_COPYLEFT_CGAL_MINKOWSKI_SUM_H
#define IGL_COPYLEFT_CGAL_MINKOWSKI_SUM_H

#include "../../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  namespace copyleft
  {
    namespace boolean
    {
      // Compute the Minkowski sum of a closed triangle mesh (V,F) and a
      // segment [s,d] in 3D.
      //
      // Inputs:
      //   VA  #VA by 3 list of mesh vertices in 3D
      //   FA  #FA by 3 list of triangle indices into VA
      //   s  segment source endpoint in 3D
      //   d  segment source endpoint in 3D
      //   resolve_overlaps  whether or not to resolve self-union. If false
      //     then result may contain self-intersections if input mesh is
      //     non-convex.
      // Outputs:
      //   W  #W by 3 list of mesh vertices in 3D
      //   G  #G by 3 list of triangle indices into W
      //   J  #G list of indices into [F;#V+F;[s d]] of birth parents
      //
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename Deriveds,
        typename Derivedd,
        typename DerivedW,
        typename DerivedG,
        typename DerivedJ>
      IGL_INLINE void minkowski_sum(
        const Eigen::PlainObjectBase<DerivedVA> & VA,
        const Eigen::PlainObjectBase<DerivedFA> & FA,
        const Eigen::PlainObjectBase<Deriveds> & s,
        const Eigen::PlainObjectBase<Derivedd> & d,
        const bool resolve_overlaps,
        Eigen::PlainObjectBase<DerivedW> & W,
        Eigen::PlainObjectBase<DerivedG> & G,
        Eigen::PlainObjectBase<DerivedJ> & J);
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename Deriveds,
        typename Derivedd,
        typename DerivedW,
        typename DerivedG,
        typename DerivedJ>
      IGL_INLINE void minkowski_sum(
        const Eigen::PlainObjectBase<DerivedVA> & VA,
        const Eigen::PlainObjectBase<DerivedFA> & FA,
        const Eigen::PlainObjectBase<Deriveds> & s,
        const Eigen::PlainObjectBase<Derivedd> & d,
        Eigen::PlainObjectBase<DerivedW> & W,
        Eigen::PlainObjectBase<DerivedG> & G,
        Eigen::PlainObjectBase<DerivedJ> & J);

    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "minkowski_sum.cpp"
#endif

#endif
