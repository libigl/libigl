//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef ISBORDERVERTEX_H
#define ISBORDERVERTEX_H

#include <Eigen/Core>
#include <string>

#include <vector>
#include "tt.h"

namespace igl 
{
  template<typename T>
  inline std::vector<bool> is_border_vertex(const T& V, const Eigen::MatrixXi& F)
  {
    Eigen::MatrixXi FF;
    igl::tt(V,F,FF);
    vector<bool> ret(V.rows());
    for(int i=0; i<ret.size();++i)
      ret[i] = false;
    
    for(int i=0; i<F.rows();++i)
      for(int j=0;j<F.cols();++j)
        if(FF(i,j) == -1)
        {
          ret[F(i,j)]       = true;
          ret[F(i,(j+1)%3)] = true;
        }
    return ret;
  }
}

#endif
