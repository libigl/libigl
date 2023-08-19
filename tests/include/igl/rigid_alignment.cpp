// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/rigid_alignment.h>
#include <igl/matlab_format.h>
#include <iostream>


TEST_CASE("rigid_alignment: identity", "[igl]")
{
  const Eigen::MatrixXd X = (Eigen::MatrixXd(10,3)<<
    0.814724,0.157613,0.655741,
    0.905792,0.970593,0.035712,
    0.126987,0.957167,0.849129,
    0.913376,0.485376,0.933993,
    0.632359,0.800280,0.678735,
    0.097540,0.141886,0.757740,
    0.278498,0.421761,0.743132,
    0.546882,0.915736,0.392227,
    0.957507,0.792207,0.655478,
    0.964889,0.959492,0.171187).finished();
  const Eigen::MatrixXd Y = X;
  const Eigen::MatrixXd N = (Y.array()-0.5).matrix().rowwise().normalized();
  Eigen::Matrix3d R;
  Eigen::RowVector3d t;
  igl::rigid_alignment(X,Y,N,R,t);
  test_common::assert_near(R,Eigen::Matrix3d::Identity(),1e-12);
  test_common::assert_near(t,Eigen::RowVector3d::Zero(),1e-12);
}
