#ifndef IGL_ORIENT_OUTWARD_AO_H
#define IGL_ORIENT_OUTWARD_AO_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Forward define
  template <
    typename PointMatrixType,
    typename FaceMatrixType,
    typename RowVector3>
  class EmbreeIntersector;
  // Orient each component (identified by C) of a mesh (V,F) using ambient occlusion 
  // such that the front side is less occluded than back side
  //
  // Inputs:
  //   V            #V by 3 list of vertex positions
  //   F            #F by 3 list of triangle indices
  //   C            #F list of components
  //   ei           EmbreeIntersector containing (V,F)
  //   num_samples  total number of rays to be shot
  // Outputs:
  //   FF  #F by 3 list of new triangle indices such that FF(~I,:) = F(~I,:) and
  //     FF(I,:) = fliplr(F(I,:)) (OK if &FF = &F)
  //   I  max(C)+1 list of whether face has been flipped
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedC, 
    typename PointMatrixType,
    typename FaceMatrixType,
    typename RowVector3,
    typename DerivedFF, 
    typename DerivedI>
  IGL_INLINE void orient_outward_ao(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedC> & C,
    const igl::EmbreeIntersector<PointMatrixType,FaceMatrixType,RowVector3> & ei,
    const int num_samples,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
  
  // EmbreeIntersector generated on the fly
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedC, 
    typename DerivedFF, 
    typename DerivedI>
  IGL_INLINE void orient_outward_ao(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedC> & C,
    const int num_samples,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
};

#ifdef IGL_HEADER_ONLY
#  include "orient_outward_ao.cpp"
#endif

#endif
