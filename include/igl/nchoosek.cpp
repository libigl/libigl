// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "nchoosek.h"

namespace igl {
  class CombinationFinder
  {
  private:
    std::vector<int> combinations;
    void add(const std::vector<int>& v,
             std::vector<std::vector<int> > &allCombs)
    {
      allCombs.push_back(v);
    }

  public:
    void doCombs(int offset,
                 int k,
                 int N,
                 std::vector<std::vector<int> > &allCombs)
    {
      if (k == 0) {
        add(combinations,allCombs);
        return;
      }
      for (int i = offset; i <= N - k; ++i) {
        combinations.push_back(i);
        doCombs(i+1, k-1, N,allCombs);
        combinations.pop_back();
      }
    }

  };


}

IGL_INLINE void igl::nchoosek(int offset,
                              int k,
                              int N,
                              std::vector<std::vector<int> > &allCombs)
{
  CombinationFinder cmbf;
  allCombs.clear();
  cmbf.doCombs(offset,k,N,allCombs);
}



#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
