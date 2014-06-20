#ifndef IGL_COMPUTE_FRAME_FIELD_BISECTORS_H
#define IGL_COMPUTE_FRAME_FIELD_BISECTORS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute face normals via vertex position list, face list
  // Inputs:
  //   V     #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F     #F by 3 eigen Matrix of face (triangle) indices
  //   B1    #F by 3 eigen Matrix of face (triangle) base vector 1
  //   B2    #F by 3 eigen Matrix of face (triangle) base vector 2
  //   PD1   #F by 3 eigen Matrix of the first per face frame field vector
  //   PD2   #F by 3 eigen Matrix of the second per face frame field vector
  // Output:
  //   BIS1  #F by 3 eigen Matrix of the first per face frame field bisector
  //   BIS2  #F by 3 eigen Matrix of the second per face frame field bisector
  //
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void compute_frame_field_bisectors(
                                                const Eigen::PlainObjectBase<DerivedV>& V,
                                                const Eigen::PlainObjectBase<DerivedF>& F,
                                                const Eigen::PlainObjectBase<DerivedV>& B1,
                                                const Eigen::PlainObjectBase<DerivedV>& B2,
                                                const Eigen::PlainObjectBase<DerivedV>& PD1,
                                                const Eigen::PlainObjectBase<DerivedV>& PD2,
                                                Eigen::PlainObjectBase<DerivedV>& BIS1,
                                                Eigen::PlainObjectBase<DerivedV>& BIS2);
  // Wrapper without given basis vectors.
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void compute_frame_field_bisectors(
                                                const Eigen::PlainObjectBase<DerivedV>& V,
                                                const Eigen::PlainObjectBase<DerivedF>& F,
                                                const Eigen::PlainObjectBase<DerivedV>& PD1,
                                                const Eigen::PlainObjectBase<DerivedV>& PD2,
                                                Eigen::PlainObjectBase<DerivedV>& BIS1,
                                                Eigen::PlainObjectBase<DerivedV>& BIS2);
}

#ifdef IGL_HEADER_ONLY
#  include "compute_frame_field_bisectors.cpp"
#endif

#endif
