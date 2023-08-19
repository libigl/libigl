// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/polygon_corners.h>
#include <igl/matrix_to_list.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("polygon_corners: quads", "[igl]")
{
  const Eigen::MatrixXi Q = (Eigen::MatrixXi(2,4)<< 0,1,2,3, 1,4,5,2).finished();
  std::vector<std::vector<int>> vQ;
  igl::matrix_to_list(Q,vQ);
  const Eigen::VectorXi Iexact = (Eigen::VectorXi(8)<<0,1,2,3,1,4,5,2).finished();
  const Eigen::VectorXi Cexact = (Eigen::VectorXi(3)<<0,4,8).finished();
  Eigen::VectorXi I,C;
  igl::polygon_corners(vQ,I,C);
  test_common::assert_eq(I,Iexact);
  test_common::assert_eq(C,Cexact);
  igl::polygon_corners( Q,I,C);
  test_common::assert_eq(I,Iexact);
  test_common::assert_eq(C,Cexact);
}
