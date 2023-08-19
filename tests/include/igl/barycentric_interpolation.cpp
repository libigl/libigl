// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/barycentric_interpolation.h>

TEST_CASE("barycentric_interpolation: two-triangles", "[igl]")
{
  const Eigen::MatrixXd V = 
    (Eigen::MatrixXd(4,3) << 
     0,0,0,
     1,0,0,
     1,1,0,
     0,1,0).finished();
  const Eigen::MatrixXi F = 
    (Eigen::MatrixXi(2,3) << 0,1,2,0,2,3).finished();
  const Eigen::VectorXi I =
   (Eigen::VectorXi(4,1)<<0,0,1,1).finished();
  const Eigen::MatrixXd B =
   (Eigen::MatrixXd(4,3)<<
    1,0,0,
    0.5,0.5,0,
    0.5,0.5,0,
    0.25,0.25,0.5
    ).finished();
  const Eigen::MatrixXd Xgt =
   (Eigen::MatrixXd(4,3)<<
    0,0,0,
    0.5,0,0,
    0.5,0.5,0,
    0.25,0.75,0
    ).finished();
  Eigen::MatrixXd X;
  igl::barycentric_interpolation(V,F,B,I,X);
  test_common::assert_eq(X,Xgt);
}
