//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_IS_BORDER_VERTEX_H
#define IGL_IS_BORDER_VERTEX_H

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  template<typename T>
  inline std::vector<bool> is_border_vertex(const T& V, const Eigen::MatrixXi& F)
}

// Implementation
#include "tt.h"

template<typename T>
inline std::vector<bool> igl::is_border_vertex(const T& V, const Eigen::MatrixXi& F)
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


#endif
