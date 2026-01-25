#include "cubic.h"
  template
    <
    typename  DerivedC,
    typename  DerivedP>
  IGL_INLINE void igl::cubic(
    const Eigen::MatrixBase<DerivedC>& C,
    const typename DerivedP::Scalar & t,
    Eigen::PlainObjectBase<DerivedP>& P)
{
  using Scalar = typename DerivedP::Scalar;
  // Evaluate cubic Bezier at parameter t
  P = 
    (Scalar(1) - t) * (Scalar(1) - t) * (Scalar(1) - t) * C.row(0)
    + Scalar(3) * (Scalar(1) - t) * (Scalar(1) - t) * t * C.row(1)
    + Scalar(3) * (Scalar(1) - t) * t * t * C.row(2)
    + t * t * t * C.row(3);
}



#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::cubic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, -1, 1, 1, -1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::Matrix<double, 1, -1, 1, 1, -1>::Scalar const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1>>&);
template void igl::cubic<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, -1, 1, true>, Eigen::Matrix<double, 1, 1, 0, 1, 1>>(Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, -1, 1, true>> const&, Eigen::Matrix<double, 1, 1, 0, 1, 1>::Scalar const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 1, 0, 1, 1>>&);
template void igl::cubic<Eigen::IndexedView<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Block<Eigen::Matrix<int, -1, -1, 0, -1, -1>, 1, -1, false>, Eigen::internal::AllRange<-1>>, Eigen::Matrix<double, 1, -1, 1, 1, -1>>(Eigen::MatrixBase<Eigen::IndexedView<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Block<Eigen::Matrix<int, -1, -1, 0, -1, -1>, 1, -1, false>, Eigen::internal::AllRange<-1>>> const&, Eigen::Matrix<double, 1, -1, 1, 1, -1>::Scalar const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1>>&);
#endif
