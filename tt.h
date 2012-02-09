//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_TT_H
#define IGL_TT_H
#include "igl_inline.h"

#include <Eigen/Core>

#include <vector>

namespace igl 
{
  // Preprocessing
  template<typename T, typename S> 
  IGL_INLINE void tt_preprocess(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& F, std::vector<std::vector<int> >& TTT);
  // Extract the face adjacencies
  template<typename S> 
  IGL_INLINE void tt_extractTT(const Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& F, std::vector<std::vector<int> >& TTT, Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& TT);
  // Extract the face adjacencies indices (needed for fast traversal)
  template<typename S> 
  IGL_INLINE void tt_extractTTi(const Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& F, std::vector<std::vector<int> >& TTT, Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& TTi);
  // Compute triangle-triangle adjacency
  template<typename T, typename S> 
  IGL_INLINE void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& F, Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& TT);
  // Compute triangle-triangle adjacency with indices
  template<typename T, typename S> 
  IGL_INLINE void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& F, Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& TT, Eigen::Matrix<S,Eigen::Dynamic, Eigen::Dynamic>& TTi);
}

#ifdef IGL_HEADER_ONLY
#  include "tt.cpp"
#endif

#endif
