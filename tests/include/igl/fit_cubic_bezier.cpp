// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/fit_cubic_bezier.h>
#include <igl/bezier.h>
#include <igl/PI.h>

TEST_CASE("fit_cubic_bezier: hemicircle", "[igl]")
{
  // Create a hemicircle
  Eigen::MatrixXd d(101,2);
  for(int i=0;i<d.rows();i++)
  {
    const double th = double(i)/double(d.rows()-1)*igl::PI;
    d(i,0) = cos(th);
    d(i,1) = sin(th);
  }
  const double error = 0.000001;
  std::vector<Eigen::MatrixXd> cubics;
  igl::fit_cubic_bezier(d,error,cubics);
  REQUIRE(cubics.size()>1);
  REQUIRE(cubics.size()<10);
  // Generate a dense sampling
  const Eigen::VectorXd T = Eigen::VectorXd::LinSpaced(1000,0,1);
  Eigen::MatrixXd X;
  igl::bezier(cubics,T,X);
  for(int j=0;j<d.rows();j++)
  {
    const double sd = 
      ((X.rowwise()-d.row(j)).rowwise().squaredNorm()).minCoeff();
    REQUIRE(sd < error);
  }
}
