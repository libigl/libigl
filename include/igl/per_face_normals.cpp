// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "per_face_normals.h"
#include <Eigen/Geometry>

#define SQRT_ONE_OVER_THREE 0.57735026918962573
template <typename DerivedV, typename DerivedF, typename DerivedZ, typename DerivedN>
IGL_INLINE void igl::per_face_normals(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const Eigen::PlainObjectBase<DerivedZ> & Z,
  Eigen::PlainObjectBase<DerivedN> & N)
{
  N.resize(F.rows(),3);
  // loop over faces
  int Frows = F.rows();
#pragma omp parallel for if (Frows>10000)
  for(int i = 0; i < Frows;i++)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = V.row(F(i,1)) - V.row(F(i,0));
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = V.row(F(i,2)) - V.row(F(i,0));
    N.row(i) = v1.cross(v2);//.normalized();
    typename DerivedV::Scalar r = N.row(i).norm();
    if(r == 0)
    {
      N.row(i) = Z;
    }else
    {
      N.row(i) /= r;
    }
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedN>
IGL_INLINE void igl::per_face_normals(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedN> & N)
{
  using namespace Eigen;
  Matrix<typename DerivedN::Scalar,3,1> Z(0,0,0);
  return per_face_normals(V,F,Z,N);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::per_face_normals<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::per_face_normals<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::per_face_normals<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&);
template void igl::per_face_normals<Eigen::Matrix<float, -1, -1, 1, -1, -1>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1>, Eigen::Matrix<float, -1, -1, 1, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 1, -1, -1> >&);
template void igl::per_face_normals<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&);
#endif
