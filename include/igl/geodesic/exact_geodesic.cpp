// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "exact_geodesic.h"
#include <geodesic_algorithm_exact.h>
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedVS,
  typename DerivedFS,
  typename DerivedVT,
  typename DerivedFT,
  typename DerivedD>
IGL_INLINE void igl::geodesic::exact_geodesic(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedVS> &VS,
  const Eigen::MatrixBase<DerivedFS> &FS,
  const Eigen::MatrixBase<DerivedVT> &VT,
  const Eigen::MatrixBase<DerivedFT> &FT,
  Eigen::PlainObjectBase<DerivedD> &D)
{
  assert(V.cols() == 3 && F.cols() == 3 && "Only support 3D triangle mesh");
  assert(VS.cols() ==1 && FS.cols() == 1 && VT.cols() == 1 && FT.cols() ==1 && "Only support one dimensional inputs");
  std::vector<typename DerivedV::Scalar> points(V.rows() * V.cols());
  std::vector<typename DerivedF::Scalar> faces(F.rows() * F.cols());
  for (int i = 0; i < points.size(); i++)
  {
    points[i] = V(i / 3, i % 3);
  }
  for (int i = 0; i < faces.size(); i++)
  {
    faces[i] = F(i / 3, i % 3);
  }

  ::geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  ::geodesic::GeodesicAlgorithmExact exact_algorithm(&mesh);

  std::vector<::geodesic::SurfacePoint> source(VS.rows() + FS.rows());
  std::vector<::geodesic::SurfacePoint> target(VT.rows() + FT.rows());
  for (int i = 0; i < VS.rows(); i++)
  {
    source[i] = (::geodesic::SurfacePoint(&mesh.vertices()[VS(i)]));
  }
  for (int i = 0; i < FS.rows(); i++)
  {
    source[i] = (::geodesic::SurfacePoint(&mesh.faces()[FS(i)]));
  }

  for (int i = 0; i < VT.rows(); i++)
  {
    target[i] = (::geodesic::SurfacePoint(&mesh.vertices()[VT(i)]));
  }
  for (int i = 0; i < FT.rows(); i++)
  {
    target[i] = (::geodesic::SurfacePoint(&mesh.faces()[FT(i)]));
  }

  exact_algorithm.propagate(source);
  std::vector<::geodesic::SurfacePoint> path;
  D.resize(target.size(), 1);
  for (int i = 0; i < target.size(); i++)
  {
    exact_algorithm.trace_back(target[i], path);
    D(i) = ::geodesic::length(path);
  }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::geodesic::exact_geodesic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>> &);
template void igl::geodesic::exact_geodesic< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(
       Eigen::MatrixBase<Eigen::Matrix<double,-1,-1,0,-1,-1>const &, 
       Eigen::MatrixBase<Eigen::Matrix<int,-1,-1,0,-1,-1> const &, 
       Eigen::MatrixBase<Eigen::Matrix<int,-1,-1,0,-1,-1> const &, 
       Eigen::MatrixBase<Eigen::Matrix<int,-1,-1,0,-1,-1> const &, 
       Eigen::MatrixBase<Eigen::Matrix<int,-1,-1,0,-1,-1> const &, 
       Eigen::MatrixBase<Eigen::Matrix<int,-1,-1,0,-1,-1> const &, 
       Eigen::PlainObjectBase<Eigen::Matrix<double,-1,-1,0,-1,-1>&);

#endif
