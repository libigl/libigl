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
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void tt_preprocess(const Eigen::PlainObjectBase<DerivedV>& V,
                                const Eigen::PlainObjectBase<DerivedF>& F,
                                std::vector<std::vector<int> >& TTT);
  // Extract the face adjacencies
  template <typename DerivedF, typename DerivedTT>
  IGL_INLINE void tt_extractTT(const Eigen::PlainObjectBase<DerivedF>& F,
                               std::vector<std::vector<int> >& TTT,
                               Eigen::PlainObjectBase<DerivedTT>& TT);
  // Extract the face adjacencies indices (needed for fast traversal)
  template <typename DerivedF, typename DerivedTT>
  IGL_INLINE void tt_extractTTi(const Eigen::PlainObjectBase<DerivedF>& F,
                                std::vector<std::vector<int> >& TTT,
                                Eigen::PlainObjectBase<DerivedTT>& TTi);
  // Compute triangle-triangle adjacency
  template <typename DerivedV, typename DerivedF, typename DerivedTT>
  IGL_INLINE void tt(const Eigen::PlainObjectBase<DerivedV>& V,
                     const Eigen::PlainObjectBase<DerivedF>& F,
                     Eigen::PlainObjectBase<DerivedTT>& TT);
  // Compute triangle-triangle adjacency with indices
  template <typename DerivedV, typename DerivedF, typename DerivedTT>
  IGL_INLINE void tt(const Eigen::PlainObjectBase<DerivedV>& V,
                     const Eigen::PlainObjectBase<DerivedF>& F,
                     Eigen::PlainObjectBase<DerivedTT>& TT,
                     Eigen::PlainObjectBase<DerivedTT>& TTi);
}

#ifdef IGL_HEADER_ONLY
#  include "tt.cpp"
#endif

#endif
