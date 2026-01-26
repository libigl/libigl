#ifndef IGL_CYCODEBASE_POINT_CUBIC_SQUARED_DISTANCE_H
#define IGL_CYCODEBASE_POINT_CUBIC_SQUARED_DISTANCE_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
  namespace cycodebase {
    template <
      typename DerivedQ,
      typename DerivedC,
      typename DerivedsqrD,
      typename DerivedS,
      typename DerivedK>
    void point_cubic_squared_distance(
      const Eigen::MatrixBase<DerivedQ>& Q,
      const Eigen::MatrixBase<DerivedC>& C,
      Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
      Eigen::PlainObjectBase<DerivedS>& S,
      Eigen::PlainObjectBase<DerivedK>& K);
    template <
      typename DerivedQ,
      typename DerivedC,
      typename DerivedD,
      typename DerivedB,
      typename DerivedsqrD,
      typename DerivedS,
      typename DerivedK>
    void point_cubic_squared_distance(
      const Eigen::MatrixBase<DerivedQ>& Q,
      const Eigen::MatrixBase<DerivedC>& C,
      const Eigen::MatrixBase<DerivedD>& D,
      const Eigen::MatrixBase<DerivedB>& B,
      Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
      Eigen::PlainObjectBase<DerivedS>& S,
      Eigen::PlainObjectBase<DerivedK>& K);
    template <
      typename Derivedq,
      typename DerivedC,
      typename DerivedD,
      typename DerivedB,
      typename Derivedk
      >
    void point_cubic_squared_distance(
      const Eigen::MatrixBase<Derivedq>& q,
      const Eigen::MatrixBase<DerivedC>& C,
      const Eigen::MatrixBase<DerivedD>& D,
      const Eigen::MatrixBase<DerivedB>& B,
      typename Derivedq::Scalar& sqrD,
      typename Derivedq::Scalar& s,
      Eigen::PlainObjectBase<Derivedk>& k);
    template <
      typename Derivedq,
      typename DerivedC,
      typename Derivedk
      >
    void point_cubic_squared_distance(
      const Eigen::MatrixBase<Derivedq>& q,
      const Eigen::MatrixBase<DerivedC>& C,
      typename Derivedq::Scalar& sqrD,
      typename Derivedq::Scalar& s,
      Eigen::PlainObjectBase<Derivedk>& k);
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "point_cubic_squared_distance.cpp"
#endif

#endif


