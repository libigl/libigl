// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Gustavo Segovia <gustavo.segovia@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "segments_to_arrows.h"

IGL_INLINE void igl::segments_to_arrows(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, Eigen::MatrixXd& A1, Eigen::MatrixXd& A2)
{
  Eigen::MatrixXd vector = (P2 - P1) * (4.0/5);
  Eigen::MatrixXd vector_perpend(vector.rows(), vector.cols());
  vector_perpend.col(0) = - vector.col(1);
  vector_perpend.col(1) =   vector.col(0);
  vector_perpend.col(2) = Eigen::VectorXd::Zero(vector.rows());

  Eigen::MatrixXd arrow_head_P1_side1 = vector + vector_perpend * 0.1 + P1;
  Eigen::MatrixXd arrow_head_P1_side2 = vector - vector_perpend * 0.1 + P1;

  A1.resize(P1.rows() * 3, P1.cols());
  A1 << P1, arrow_head_P1_side1, arrow_head_P1_side2;
  A2.resize(P2.rows() * 3, P2.cols());
  A2 << P2, P2, P2;
}

