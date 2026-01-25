#ifndef IGL_CYCODEBASE_POINT_SPLINE_SQUARED_DISTANCE_H
#define IGL_CYCODEBASE_POINT_SPLINE_SQUARED_DISTANCE_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
  namespace cycodebase {
    template <
      typename DerivedQ,
      typename DerivedP,
      typename DerivedC,
      typename DerivedsqrD,
      typename DerivedI,
      typename DerivedS,
      typename DerivedK>
    void point_spline_squared_distance(
      const Eigen::MatrixBase<DerivedQ>& Q,
      const Eigen::MatrixBase<DerivedP>& P,
      const Eigen::MatrixBase<DerivedC>& C,
      Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
      Eigen::PlainObjectBase<DerivedI>& I,
      Eigen::PlainObjectBase<DerivedS>& S,
      Eigen::PlainObjectBase<DerivedK>& K);
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "point_spline_squared_distance.cpp"
#endif

#endif



