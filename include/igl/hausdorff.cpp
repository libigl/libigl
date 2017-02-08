// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "hausdorff.h"
#include "point_mesh_squared_distance.h"

template <
  typename DerivedVA, 
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename Scalar>
IGL_INLINE void igl::hausdorff(
  const Eigen::PlainObjectBase<DerivedVA> & VA, 
  const Eigen::PlainObjectBase<DerivedFA> & FA,
  const Eigen::PlainObjectBase<DerivedVB> & VB, 
  const Eigen::PlainObjectBase<DerivedFB> & FB,
  Scalar & d)
{
  using namespace Eigen;
  assert(VA.cols() == 3 && "VA should contain 3d points");
  assert(FA.cols() == 3 && "FA should contain triangles");
  assert(VB.cols() == 3 && "VB should contain 3d points");
  assert(FB.cols() == 3 && "FB should contain triangles");
  Matrix<Scalar,Dynamic,1> sqr_DBA,sqr_DAB;
  Matrix<typename DerivedVA::Index,Dynamic,1> I;
  Matrix<typename DerivedVA::Scalar,Dynamic,3> C;
  point_mesh_squared_distance(VB,VA,FA,sqr_DBA,I,C);
  point_mesh_squared_distance(VA,VB,FB,sqr_DAB,I,C);
  const Scalar dba = sqr_DBA.maxCoeff();
  const Scalar dab = sqr_DAB.maxCoeff();
  d = sqrt(std::max(dba,dab));
}

#ifdef IGL_STATIC_LIBRARY
template void igl::hausdorff<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, double&);
#endif
