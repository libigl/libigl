#include "edge_lengths.h"

template <typename DerivedV, typename DerivedF, typename DerivedL>
IGL_INLINE void igl::edge_lengths(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedL>& L)
{
  assert(F.cols() == 3);
  L.resize(F.rows(),3);
  // loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    L(i,0) = sqrt((V.row(F(i,1))-V.row(F(i,2))).array().pow(2).sum());
    L(i,1) = sqrt((V.row(F(i,2))-V.row(F(i,0))).array().pow(2).sum());
    L(i,2) = sqrt((V.row(F(i,0))-V.row(F(i,1))).array().pow(2).sum());
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::edge_lengths<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
