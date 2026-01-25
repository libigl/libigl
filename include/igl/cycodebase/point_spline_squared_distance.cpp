#include "point_spline_squared_distance.h"
#include "box_cubic.h"

template <
  typename DerivedQ,
  typename DerivedP,
  typename DerivedC,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedS,
  typename DerivedK>
void igl::cycodebase::point_spline_squared_distance(
  const Eigen::MatrixBase<DerivedQ>& Q,
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
  Eigen::PlainObjectBase<DerivedI>& I,
  Eigen::PlainObjectBase<DerivedS>& S,
  Eigen::PlainObjectBase<DerivedK>& K)
{
  using Scalar = typename DerivedQ::Scalar;
  Eigen::Matrix<Scalar,DerivedC::RowsAtCompileTime,DerivedP::ColsAtCompileTime> B1,B2;
  box_cubic(P,C,B1,B2);
  //point_cubic_squared_distance(
}
