#include "histc.h"
#include <cassert>

template <typename DerivedX, typename DerivedE, typename DerivedN, typename DerivedB>
IGL_INLINE void igl::histc(
  const Eigen::PlainObjectBase<DerivedX > & X,
  const Eigen::PlainObjectBase<DerivedE > & E,
  Eigen::PlainObjectBase<DerivedN > & N,
  Eigen::PlainObjectBase<DerivedB > & B)
{
  const int m = X.size();
  const int n = E.size();
  assert( 
    (E.topLeftCorner(E.size()-1) - 
      E.bottomRightCorner(E.size()-1)).maxCoeff() >= 0 && 
    "E should be monotonically increasing");
  N.resize(n,1);
  B.resize(m,1);
  N.setConstant(0);
#pragma omp parallel for
  for(int j = 0;j<m;j++)
  {
    const double x = X(j);
    // Boring one-offs
    if(x < E(0) || x > E(E.size()-1))
    {
      B(j) = -1;
      continue;
    }
    // Find x in E
    int l = 0;
    int h = E.size()-1;
    int k = l;
    while((h-l)>1)
    {
      assert(x >= E(l));
      assert(x <= E(h));
      k = (h+l)/2;
      if(x < E(k))
      {
        h = k;
      }else
      {
        l = k;
      }
    }
    if(x == E(h))
    {
      k = h;
    }else
    {
      k = l;
    }
    B(j) = k;
    // Not sure if this pragma is needed
#pragma omp atomic
    N(k)++;
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
