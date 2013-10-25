#include "randperm.h"
#include "colon.h"
#include <algorithm> 

template <typename DerivedI>
IGL_INLINE void igl::randperm(
  const int n,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  Eigen::VectorXi II;
  igl::colon(0,1,n-1,II);
  I = II;
  std::random_shuffle(I.data(),I.data()+n);
}

template <typename DerivedI>
IGL_INLINE Eigen::PlainObjectBase<DerivedI> igl::randperm( const int n)
{
  Eigen::PlainObjectBase<DerivedI> I;
  randperm(n,I);
  return I;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
