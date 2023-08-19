// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/blkdiag.h>

TEST_CASE("blkdiag: 3-matrices", "[igl]")
{
  const Eigen::MatrixXd Ygt= 
    (Eigen::MatrixXd(3+0+2,3+2+5)<<
      1,2,3,0,0,0,0,0,0,0,
      4,5,6,0,0,0,0,0,0,0,
      7,8,9,0,0,0,0,0,0,0,
      0,0,0,0,0,0,1,2,3,4,
      0,0,0,0,0,5,6,7,8,9).finished();
  Eigen::MatrixXd A(3,3);
  A<<1,2,3,4,5,6,7,8,9;
  // This is the correct behavior for a 0×n matrix
  Eigen::MatrixXd B(0,2);
  Eigen::MatrixXd C(2,5);
  C<<0,1,2,3,4,5,6,7,8,9;
  Eigen::MatrixXd Y;
  igl::blkdiag({A,B,C},Y);
  test_common::assert_eq(Y,Ygt);
  {
    Eigen::SparseMatrix<double> sA = A.sparseView();
    Eigen::SparseMatrix<double> sB = B.sparseView();
    Eigen::SparseMatrix<double> sC = C.sparseView();
    Eigen::SparseMatrix<double> sY;
    igl::blkdiag({sA,sB,sC},sY);
    Y = sY;
    test_common::assert_eq(Y,Ygt);
  }
}

