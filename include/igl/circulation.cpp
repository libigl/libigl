// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "circulation.h"
#include "list_to_matrix.h"
#include <cassert>

IGL_INLINE std::vector<int> igl::circulation(
  const int e,
  const bool ccw,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI)
{
  // prepare output
  std::vector<int> N;
  N.reserve(6);
  const int m = EMAP.size()/3;
  assert(m*3 == EMAP.size());
  const auto & step = [&](
    const int e, 
    const int ff,
    int & ne, 
    int & nf)
  {
    assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
    //const int fside = EF(e,1)==ff?1:0;
    const int nside = EF(e,0)==ff?1:0;
    const int nv = EI(e,nside);
    // get next face
    nf = EF(e,nside);
    // get next edge 
    const int dir = ccw?-1:1;
    ne = EMAP(nf+m*((nv+dir+3)%3));
  };
  // Always start with first face (ccw in step will be sure to turn right
  // direction)
  const int f0 = EF(e,0);
  int fi = f0;
  int ei = e;
  while(true)
  {
    step(ei,fi,ei,fi);
    N.push_back(fi);
    // back to start?
    if(fi == f0)
    {
      assert(ei == e);
      break;
    }
  }
  return N;
}

IGL_INLINE void igl::circulation(
  const int e,
  const bool ccw,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  Eigen::VectorXi & vN)
{
  std::vector<int> N = circulation(e,ccw,EMAP,EF,EI);
  igl::list_to_matrix(N,vN);
}

IGL_INLINE void igl::circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  /*std::vector<int> & Ne,*/
  std::vector<int> & Nv,
  std::vector<int> & Nf)
{
  //
  // for e --> (bf) and ccw=true
  //
  //     c---d
  //    / \ / \
  //   a---b-e-f
  //    \ / \ /
  //     g---h
  //
  //  // (might start with {bhf} depending on edge)
  //  Ne = […] -> [fd db dc cb ca ab ag gb gh hb hf fb]
  //              {upto cylic order}
  //  Nf = […] -> [{bfd}, {bdc}, {bca}, {bag}, {bgh}, {bhf}]
  //  Nv = [d c a g h f]
  //
  // prepare output
  //Ne.clear();Ne.reserve(2*10);
  Nv.clear();Nv.reserve(10);
  Nf.clear();Nf.reserve(10);
  const int m = EMAP.size()/3;
  assert(m*3 == EMAP.size());
  const auto & step = [&](
    const int e, 
    const int ff,
    int & ne, 
    //int & re,
    int & rv,
    int & nf)
  {
    assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
    //const int fside = EF(e,1)==ff?1:0;
    const int nside = EF(e,0)==ff?1:0;
    const int nv = EI(e,nside);
    // get next face
    nf = EF(e,nside);
    // get next edge 
    const int dir = ccw?-1:1;
    rv = F(nf,nv);
    ne = EMAP(nf+m*((nv+dir+3)%3));
    //re = EMAP(nf+m*((nv+2*dir+3)%3));
  };
  // Always start with first face (ccw in step will be sure to turn right
  // direction)
  const int f0 = EF(e,0);
  int fi = f0;
  int ei = e;
  while(true)
  {
    int re,rv;
    step(ei,fi,ei/*,re*/,rv,fi);
    Nf.push_back(fi);
    //Ne.push_back(re);
    //Ne.push_back(ei);
    Nv.push_back(rv);
    // back to start?
    if(fi == f0)
    {
      assert(ei == e);
      break;
    }
  }
}
