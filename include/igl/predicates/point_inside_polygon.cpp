// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Hanxiao Shen <hanxiao@cs.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "point_inside_polygon.h"

template <typename DerivedP, typename DerivedQ>
IGL_INLINE bool igl::predicates::point_inside_polygon(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedQ>& q
){
  for(int i=0;i<P.rows();i++){
    int i_1 = (i+1) % P.rows();
    Eigen::RowVector2d a = P.row(i);
    Eigen::RowVector2d b = P.row(i_1);
    auto r = igl::predicates::orient2d(a,b,q);
    if(r == igl::predicates::Orientation::COLLINEAR || 
       r == igl::predicates::Orientation::NEGATIVE)
      return false;
  }
  return true;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::predicates::point_inside_polygon<Eigen::Matrix<double, -1, 2, 0, -1, 2>, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&);
#endif