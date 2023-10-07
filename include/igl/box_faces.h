#ifndef IGL_BOX_FACES_H
#define IGL_BOX_FACES_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace igl
{
  template <typename DerivedV, typename DerivedQ>
  IGL_INLINE void box_faces(
    const Eigen::AlignedBox<typename DerivedV::Scalar,3> & box,
    const typename DerivedV::Scalar shrink,
    Eigen::PlainObjectBase<DerivedV> & P,
    Eigen::PlainObjectBase<DerivedQ> & Q);
  // Forward declaration
  template <typename DerivedV, int DIM> class AABB;
  template <
    typename DerivedV,
    typename DerivedP,
    typename DerivedQ,
    typename DerivedD >
  IGL_INLINE void box_faces(
    const igl::AABB<DerivedV,3> & tree,
    Eigen::PlainObjectBase<DerivedP> & P,
    Eigen::PlainObjectBase<DerivedQ> & Q,
    Eigen::PlainObjectBase<DerivedD> & D);
}

#ifndef IGL_STATIC_LIBRARY
#  include "box_faces.cpp"
#endif

#endif
