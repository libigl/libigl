// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/accumarray.h>

TEST_CASE("accumarray: matlab_help", "[igl]")
{
  const Eigen::VectorXd V = 
    (Eigen::VectorXd(5) << 101,102,103,104,105).finished();
  const Eigen::VectorXi S = 
    (Eigen::VectorXi(5) << 0,2,3,2,3).finished();
  Eigen::VectorXd A;
  igl::accumarray(S,V,A);
  const Eigen::VectorXd Agt  = 
    (Eigen::VectorXd(4) << 101,0,206,208).finished();
  test_common::assert_eq(A,Agt);
}

TEST_CASE("accumarray: scalar", "[igl]")
{
  const auto n = 
    (Eigen::VectorXi(5) << 1,2,2,4,5).finished();
  Eigen::VectorXi C;
  igl::accumarray(n,1,C);
  const auto Cgt = 
    (Eigen::VectorXi(6) << 0,1,2,0,1,1).finished();
  test_common::assert_eq(C,Cgt);
}
