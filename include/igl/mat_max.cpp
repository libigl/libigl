// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mat_max.h"

template <typename T>
IGL_INLINE void igl::mat_max(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
  const int dim,
  Eigen::Matrix<T,Eigen::Dynamic,1> & Y,
  Eigen::Matrix<int,Eigen::Dynamic,1> & I)
{
  assert(dim==1||dim==2);

  // output size
  int n = (dim==1?X.cols():X.rows());
  // resize output
  Y.resize(n);
  I.resize(n);

  // loop over dimension opposite of dim
  for(int j = 0;j<n;j++)
  {
    typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Index PHONY;
    typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Index i;
    T m;
    if(dim==1)
    {
      m = X.col(j).maxCoeff(&i,&PHONY);
    }else
    {
      m = X.row(j).maxCoeff(&PHONY,&i);
    }
    Y(j) = m;
    I(j) = i;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::mat_max<double>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, int, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);
#endif
