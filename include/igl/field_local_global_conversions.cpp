// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/field_local_global_conversions.h>

 template <typename DerivedG, typename DerivedL, typename DerivedB>
IGL_INLINE void igl::global2local(
const Eigen::PlainObjectBase<DerivedB>& B1,
const Eigen::PlainObjectBase<DerivedB>& B2,
const Eigen::PlainObjectBase<DerivedG>& global,
Eigen::PlainObjectBase<DerivedL>& local)
{

  int n = global.cols()/3;
  //get 2D vectors from 3D vectors
  local.resize(global.rows(), n*2);
  for (int i =0; i<global.rows(); ++i)
  {
    for (int k =0; k<n; k++)
    {
      const Eigen::Matrix<typename DerivedG::Scalar,1,3> &g = global.block(i, k*3,1,3);
      local.block(i, k*2,1,2) << g.dot(B1.row(i).template cast<typename DerivedG::Scalar>()), g.dot(B2.row(i).template cast<typename DerivedG::Scalar>()) ;
    }
  }
}

template <typename DerivedG, typename DerivedL, typename DerivedB>
IGL_INLINE void igl::local2global(
const Eigen::PlainObjectBase<DerivedB>& B1,
const Eigen::PlainObjectBase<DerivedB>& B2,
const Eigen::PlainObjectBase<DerivedL>& local,
Eigen::PlainObjectBase<DerivedG>& global)
{
  int n = local.cols()/2;
  //get 3D vectors from 2D vectors
  global.resize(local.rows(), n*3);
  for (int i =0; i<global.rows(); ++i)
  {
    for (int k =0; k<n; k++)
    {
      const Eigen::Matrix<typename DerivedL::Scalar,1,2> &l = local.block(i, k*2,1,2);
      global.block(i, k*3, 1,3) = l[0]*B1.row(i).template cast<typename DerivedG::Scalar>() + l[1]*B2.row(i).template cast<typename DerivedG::Scalar>() ;
    }
  }

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::global2local<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::local2global<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
