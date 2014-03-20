// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "forward_kinematics.h"
#include <functional>

IGL_INLINE void igl::forward_kinematics(
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  const Eigen::VectorXi & P,
  const std::vector<
    Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
  std::vector<
    Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & vQ,
  std::vector<Eigen::Vector3d> & vT)
{
  using namespace std;
  using namespace Eigen;
  const int m = BE.rows(); 
  assert(m == P.rows());
  assert(m == (int)dQ.size());
  vector<bool> computed(m,false);
  vQ.resize(m);
  vT.resize(m);
  function<void (int) > fk_helper = [&] (int b)
  {
    if(!computed[b])
    {
      if(P(b) < 0)
      {
        // base case for roots
        vQ[b] = dQ[b];
        const Vector3d r = C.row(BE(b,0)).transpose();
        vT[b] = r-dQ[b]*r;
      }else
      {
        // Otherwise first compute parent's
        const int p = P(b);
        fk_helper(p);
        vQ[b] = vQ[p] * dQ[b];
        const Vector3d r = C.row(BE(b,0)).transpose();
        vT[b] = vT[p] - vQ[b]*r + vQ[p]*r;
      }
      computed[b] = true;
    }
  };
  for(int b = 0;b<m;b++)
  {
    fk_helper(b);
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
#endif
