// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/barycentric_interpolation.h>
#include <igl/readOBJ.h>
#include <igl/blue_noise.h>
#include <igl/knn.h>
#include <igl/octree.h>
#include <igl/slice.h>

TEST_CASE("blue_noise: decimated-knight", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(test_common::data_path("decimated-knight.obj"),V,F);
  const double r = 0.01;
  Eigen::MatrixXd B,P;
  {
    Eigen::VectorXi I;
    igl::blue_noise(V,F,r,B,I,P);
  }
  // There should be ~4000 samples on this model
  REQUIRE(P.rows() > 3000);
  std::vector<std::vector<int> > point_indices;
  Eigen::MatrixXi CH;
  Eigen::MatrixXd CN;
  Eigen::VectorXd W;
  igl::octree(P,point_indices,CH,CN,W);
  Eigen::MatrixXi I;
  igl::knn(P,2,point_indices,CH,CN,W,I);
  Eigen::MatrixXd P2;
  igl::slice(P,I.col(1).eval(),1,P2);
  Eigen::VectorXd D = (P-P2).rowwise().norm();
  REQUIRE(D.minCoeff() > r);
}
