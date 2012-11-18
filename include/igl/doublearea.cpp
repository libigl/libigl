#include "doublearea.h"
#include "edge_lengths.h"
#include <cassert>

template <typename DerivedV, typename DerivedF, typename DeriveddblA>
IGL_INLINE void igl::doublearea( 
  const Eigen::PlainObjectBase<DerivedV> & V, 
  const Eigen::PlainObjectBase<DerivedF> & F, 
  Eigen::PlainObjectBase<DeriveddblA> & dblA)
{
  // Only support triangles
  assert(F.cols() == 3);
  // Compute edge lengths
  Eigen::PlainObjectBase<DerivedV> l;
  edge_lengths(V,F,l);
  return doublearea(l,dblA);
}

template <typename Derivedl, typename DeriveddblA>
IGL_INLINE void igl::doublearea( 
  const Eigen::PlainObjectBase<Derivedl> & l, 
  Eigen::PlainObjectBase<DeriveddblA> & dblA)
{
  // Only support triangles
  assert(l.cols() == 3);
  // Number of triangles
  const int m = l.rows();
  // semiperimeters
  Eigen::PlainObjectBase<Derivedl> s = l.rowwise().sum()*0.5;
  assert(s.rows() == m);
  // resize output
  dblA.resize(l.rows(),1);
  // Heron's formula for area
  for(int i = 0;i<m;i++)
  {
    dblA(i) = 2.0*sqrt(s(i)*(s(i)-l(i,0))*(s(i)-l(i,1))*(s(i)-l(i,2)));
  }
}

#ifndef IGL_HEADER_ONLY
template void igl::doublearea<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::doublearea<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
