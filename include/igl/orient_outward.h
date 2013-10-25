#ifndef IGL_ORIENT_OUTWARD_H
#define IGL_ORIENT_OUTWARD_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Orient each component (identified by C) of a mesh (V,F) so the normals on
  // average point away from the patch's centroid.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  //   C  #F list of components
  // Outputs:
  //   FF  #F by 3 list of new triangle indices such that FF(~I,:) = F(~I,:) and
  //     FF(I,:) = fliplr(F(I,:)) (OK if &FF = &F)
  //   I  max(C)+1 list of whether face has been flipped
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedC, 
    typename DerivedFF, 
    typename DerivedI>
  IGL_INLINE void orient_outward(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
};

#ifdef IGL_HEADER_ONLY
#  include "orient_outward.cpp"
#endif

#endif
