// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/iterative_closest_point.h>


TEST_CASE("iterative_closest_point: identity","[igl]" "[slow]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd VX,VY;
    Eigen::MatrixXi FX,FY;
    // Load example mesh: GetParam() will be name of mesh file
    igl::read_triangle_mesh(test_common::data_path(param), VY, FY);
    VX = VY;
    FX = FY;
    // Single iteration should find a identity
    srand(0);
    Eigen::Matrix3d R;
    Eigen::RowVector3d t;
    igl::iterative_closest_point(VX,FX,VY,FY,1000,1,R,t);
    test_common::assert_near(R,Eigen::Matrix3d::Identity(),1e-12);
    test_common::assert_near(t,Eigen::RowVector3d::Zero(),1e-12);
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}
