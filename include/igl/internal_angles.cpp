// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// Copyright (C) 2015 Daniele Panozzo <daniele.panozzo@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "internal_angles.h"
#include "edge_lengths.h"
#include "get_seconds.h"

template <typename DerivedV, typename DerivedF, typename DerivedK>
IGL_INLINE void igl::internal_angles(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  using namespace Eigen;
  using namespace std;
  if(F.cols() == 3)
  {
    // Edge lengths
    Matrix<
      typename DerivedV::Scalar,
      DerivedF::RowsAtCompileTime,
      DerivedF::ColsAtCompileTime> L;
    edge_lengths(V,F,L);

    assert(F.cols() == 3 && "F should contain triangles");
    internal_angles(L,K);
  }else
  {
    assert(V.cols() == 3 && "If F contains non-triangle facets, V must be 3D");
    K.resize(F.rows(),F.cols());
    auto corner = [](
      const Eigen::PlainObjectBase<DerivedV>& x, 
      const Eigen::PlainObjectBase<DerivedV>& y, 
      const Eigen::PlainObjectBase<DerivedV>& z)
    {
      Eigen::RowVector3d v1 = (x-y).normalized();
      Eigen::RowVector3d v2 = (z-y).normalized();
      // http://stackoverflow.com/questions/10133957/signed-angle-between-two-vectors-without-a-reference-plane
      double s = v1.cross(v2).norm();
      double c = v1.dot(v2);
      return atan2(s, c);
    };
    for(unsigned i=0; i<F.rows(); ++i)
    {
      for(unsigned j=0; j<F.cols(); ++j)
      {
        K(i,j) = corner(
            V.row(F(i,int(j-1+F.cols())%F.cols())),
            V.row(F(i,j)),
            V.row(F(i,(j+1+F.cols())%F.cols()))
            );
      }
    }
  }
}

template <typename DerivedL, typename DerivedK>
IGL_INLINE void igl::internal_angles(
  const Eigen::PlainObjectBase<DerivedL>& L,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  assert(L.cols() == 3 && "Edge-lengths should come from triangles");
  const size_t m = L.rows();
  K.resize(m,3);
  //for(int d = 0;d<3;d++)
  //{
  //  const auto & s1 = L.col(d).array();
  //  const auto & s2 = L.col((d+1)%3).array();
  //  const auto & s3 = L.col((d+2)%3).array();
  //  K.col(d) = ((s3.square() + s2.square() - s1.square())/(2.*s3*s2)).acos();
  //}
  // Minimum number of iterms per openmp thread
  #ifndef IGL_OMP_MIN_VALUE
  #  define IGL_OMP_MIN_VALUE 1000
  #endif
  #pragma omp parallel for if (m>IGL_OMP_MIN_VALUE)
  for(long long f = 0;f<m;f++)
  {
    for(size_t d = 0;d<3;d++)
    {
      const auto & s1 = L(f,d);
      const auto & s2 = L(f,(d+1)%3);
      const auto & s3 = L(f,(d+2)%3);
      K(f,d) = acos((s3*s3 + s2*s2 - s1*s1)/(2.*s3*s2));
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::internal_angles<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::internal_angles<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::internal_angles<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
