#ifndef IGL_AMBIENT_OCCLUSION_H
#define IGL_AMBIENT_OCCLUSION_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
namespace igl
{
  // Forward define
  template <
    typename Scalar,
    typename Index>
  class EmbreeIntersector;
  // Compute ambient occlusion per given point
  //
  // Inputs:
  //    ei  EmbreeIntersector containing (V,F)
  //    P  #P by 3 list of origin points
  //    N  #P by 3 list of origin normals
  // Outputs:
  //    S  #P list of ambient occlusion values between 1 (fully occluded) and 0
  //      (not occluded)
  //
  template <
    typename Scalar,
    typename Index,
    typename DerivedP,
    typename DerivedN,
    typename DerivedS >
  void ambient_occlusion(
    const igl::EmbreeIntersector<Scalar,Index> & ei,
    const Eigen::PlainObjectBase<DerivedP> & P,
    const Eigen::PlainObjectBase<DerivedN> & N,
    const int num_samples,
    Eigen::PlainObjectBase<DerivedS> & S);
  // Wrapper which builds new EmbreeIntersector for (V,F). That's expensive so
  // avoid this if repeatedly calling.
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename DerivedN,
    typename DerivedS >
  void ambient_occlusion(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedP> & P,
    const Eigen::PlainObjectBase<DerivedN> & N,
    const int num_samples,
    Eigen::PlainObjectBase<DerivedS> & S);
};
#ifdef IGL_HEADER_ONLY
#  include "ambient_occlusion.cpp"
#endif

#endif
