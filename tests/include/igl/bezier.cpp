// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/bezier.h>

TEST_CASE("bezier: ease", "[igl]")
{
  // Ease curve
  const Eigen::MatrixXd C = 
    (Eigen::MatrixXd(4,2)<<0,0,0.5,0,0.5,1,1,1).finished();
  const Eigen::VectorXd T = 
    (Eigen::VectorXd(11,1)<<0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1).finished();
  Eigen::MatrixXd P;
  igl::bezier(C,T,P);
  const Eigen::MatrixXd Pexact = (Eigen::MatrixXd(11,2)<<
    0.00000000000000000,0.00000000000000000,
    0.13600000000000001,0.02800000000000001,
    0.24800000000000008,0.10400000000000004,
    0.34199999999999997,0.21599999999999997,
    0.42400000000000004,0.35200000000000004,
    0.50000000000000000,0.50000000000000000,
    0.57600000000000007,0.64800000000000002,
    0.65799999999999992,0.78399999999999992,
    0.75200000000000011,0.89600000000000013,
    0.86400000000000010,0.97200000000000009,
    1.00000000000000000,1.00000000000000000).finished();
  test_common::assert_near(P,Pexact,1e-15);
}
