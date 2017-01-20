// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "combine.h"
#include <cassert>
template <
  typename DerivedVV, 
  typename DerivedFF, 
  typename DerivedV, 
  typename DerivedF>
IGL_INLINE void igl::combine(
  const std::vector<DerivedVV> & VV,
  const std::vector<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F)
{
  assert(VV.size() == FF.size() && 
    "Lists of verex lists and face lists should be same size");
  // Dimension of vertex positions
  const int dim = VV.size() > 0 ? VV[0].cols() : 0;
  // Simplex/element size
  const int ss = FF.size() > 0 ? FF[0].cols() : 0;
  int n = 0;
  int m = 0;
  for(int i = 0;i<VV.size();i++)
  {
    const auto & Vi = VV[i];
    const auto & Fi = FF[i];
    n+=Vi.rows();
    assert(dim == Vi.cols() && "All vertex lists should have same #columns");
    m+=Fi.rows();
    assert(ss == Fi.cols() && "All face lists should have same #columns");
  }
  V.resize(n,dim);
  F.resize(m,ss);
  {
    int kv = 0;
    int kf = 0;
    for(int i = 0;i<VV.size();i++)
    {
      const auto & Vi = VV[i];
      const int ni = Vi.rows();
      const auto & Fi = FF[i];
      const int mi = Fi.rows();
      F.block(kf,0,mi,ss) = Fi.array()+kv;
      kf+=mi;
      V.block(kv,0,ni,dim) = Vi;
      kv+=ni;
    }
    assert(kv == V.rows());
    assert(kf == F.rows());
  }
}
