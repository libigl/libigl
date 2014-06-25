// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SORTROWS_H
#define IGL_SORTROWS_H
#include "igl_inline.h"

#include <vector>
#include <Eigen/Core>
namespace igl
{
  // Act like matlab's [Y,I] = sortrows(X)
  //
  // Templates:
  //   DerivedX derived scalar type, e.g. MatrixXi or MatrixXd
  //   DerivedIX derived integer type, e.g. MatrixXi
  // Inputs:
  //   X  m by n matrix whose entries are to be sorted
  //   ascending  sort ascending (true, matlab default) or descending (false)
  // Outputs:
  //   Y  m by n matrix whose entries are sorted
  //   IX  m by n matrix of indices so that if dim = 1, then in matlab notation
  //     for j = 1:n, Y(:,j) = X(I(:,j),j); end
  template <typename DerivedX, typename DerivedIX>
  IGL_INLINE void sortrows(
    const Eigen::PlainObjectBase<DerivedX>& X,
    const bool ascending,
    Eigen::PlainObjectBase<DerivedX>& Y,
    Eigen::PlainObjectBase<DerivedIX>& IX);
}

#ifndef IGL_STATIC_LIBRARY
#  include "sortrows.cpp"
#endif

#endif

