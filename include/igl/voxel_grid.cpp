// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "voxel_grid.h"
#include "grid.h"

IGL_INLINE void igl::voxel_grid(
  const Eigen::AlignedBox3d & box, 
  const int in_s,
  const int pad_count,
  Eigen::MatrixXd & GV,
  Eigen::RowVector3i & side)
{
  using namespace Eigen;
  using namespace std;
  MatrixXd::Index si = -1;
  box.diagonal().maxCoeff(&si);
  //MatrixXd::Index si = 0;
  //assert(si>=0);
  const double s_len = box.diagonal()(si);
  assert(in_s>(pad_count*2+1) && "s should be > 2*pad_count+1");
  const double s = in_s - 2*pad_count;
  side(si) = s;
  for(int i = 0;i<3;i++)
  {
    if(i!=si)
    {
      side(i) = std::ceil(s * (box.max()(i)-box.min()(i))/s_len);
    }
  }
  side.array() += 2*pad_count;
  grid(side,GV);
  // A *    p/s  + B = min
  // A * (1-p/s) + B = max
  // B = min - A * p/s
  // A * (1-p/s) + min - A * p/s = max
  // A * (1-p/s) - A * p/s = max-min
  // A * (1-2p/s) = max-min
  // A  = (max-min)/(1-2p/s)
  const Array<double,3,1> ps= 
    (double)(pad_count)/(side.transpose().cast<double>().array()-1.);
  const Array<double,3,1> A = box.diagonal().array()/(1.0-2.*ps);
  //// This would result in an "anamorphic", but perfectly fit grid:
  //const Array<double,3,1> B = box.min().array() - A.array()*ps;
  //GV.array().rowwise() *= A.transpose();
  //GV.array().rowwise() += B.transpose();
  // Instead scale by largest factor and move to match center
  Array<double,3,1>::Index ai = -1;
  double a = A.maxCoeff(&ai);
  const Array<double,1,3> ratio = 
    a*(side.cast<double>().array()-1.0)/(double)(side(ai)-1.0);
  GV.array().rowwise() *= ratio;
  GV.rowwise() += (box.center().transpose()-GV.colwise().mean()).eval();
}
