// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "triangulate.h"

template <typename Index, typename DerivedF>
IGL_INLINE void igl::triangulate(
  const std::vector<std::vector<Index> > & vF,
  Eigen::PlainObjectBase<DerivedF>& F)
{
  using namespace std;
  using namespace Eigen;
  int m = 0;
  // estimate of size
  for(typename vector<vector<Index > >::const_iterator fit = vF.begin();
    fit!=vF.end();
    fit++)
  {
    if(fit->size() >= 3)
    {
      m += fit->size() - 2;
    }
  }
  // Resize output
  F.resize(m,3);
  {
    int k = 0;
    for(typename vector<vector<Index > >::const_iterator fit = vF.begin();
      fit!=vF.end();
      fit++)
    {
      if(fit->size() >= 3)
      {
        typename vector<Index >::const_iterator cit = fit->begin();
        cit++;
        typename vector<Index >::const_iterator pit = cit++;
        for(;
          cit!=fit->end();
          cit++,pit++)
        {
          F(k,0) = *(fit->begin());
          F(k,1) = *pit;
          F(k,2) = *cit;
          k++;
        }
      }
    }
    assert(k==m);
  }

}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::triangulate<int, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
