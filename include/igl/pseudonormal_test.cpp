// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "pseudonormal_test.h"
#include "AABB.h"
#include <cassert>

IGL_INLINE void igl::pseudonormal_test(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  const Eigen::RowVector3d & q,
  const int f,
  const Eigen::RowVector3d & c,
  double & s,
  Eigen::RowVector3d & n)
{
  using namespace Eigen;
  const auto & qc = q-c;
  RowVector3d b;
  AABB<Eigen::MatrixXd,3>::barycentric_coordinates(
    c,V.row(F(f,0)),V.row(F(f,1)),V.row(F(f,2)),b);
  // Determine which normal to use
  const double epsilon = 1e-12;
  const int type = (b.array()<=epsilon).cast<int>().sum();
  switch(type)
  {
    case 2:
      // Find vertex
      for(int x = 0;x<3;x++)
      {
        if(b(x)>epsilon)
        {
          n = VN.row(F(f,x));
          break;
        }
      }
      break;
    case 1:
      // Find edge
      for(int x = 0;x<3;x++)
      {
        if(b(x)<=epsilon)
        {
          n = EN.row(EMAP(F.rows()*x+f));
          break;
        }
      }
      break;
    default:
      assert(false && "all barycentric coords zero.");
    case 0:
      n = FN.row(f);
      break;
  }
  s = (qc.dot(n) >= 0 ? 1. : -1.);
}
