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
  template<typename T> 
  IGL_INLINE void tt_preprocess(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT);
  // Extract the face adjacencies
  IGL_INLINE void tt_extractTT(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT, Eigen::MatrixXi& TT);
  // Extract the face adjacencies indices (needed for fast traversal)
  IGL_INLINE void tt_extractTTi(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT, Eigen::MatrixXi& TTi);
  // Compute triangle-triangle adjacency
  template<typename T> 
  IGL_INLINE void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& TT);
  // Compute triangle-triangle adjacency with indices
  template<typename T> 
  IGL_INLINE void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& TT, Eigen::MatrixXi& TTi);
}

#ifdef IGL_HEADER_ONLY
#  include "tt.cpp"
#endif

#endif
