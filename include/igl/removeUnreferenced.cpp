// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "removeUnreferenced.h"

template <typename Scalar, typename Index>
IGL_INLINE void igl::removeUnreferenced(
                                        const Eigen::PlainObjectBase<Scalar> &V,
                                        const Eigen::PlainObjectBase<Index> &F,
                                        Eigen::PlainObjectBase<Scalar> &NV,
                                        Eigen::PlainObjectBase<Index> &NF,
                                        Eigen::PlainObjectBase<Index> &I
                                        )
{

  // Mark referenced vertices
  Eigen::MatrixXi mark = Eigen::MatrixXi::Zero(V.rows(),1);
  
  for(int i=0; i<F.rows(); ++i)
  {
    for(int j=0; j<F.cols(); ++j)
    {
      if (F(i,j) != -1)
        mark(F(i,j)) = 1;
    }
  }
  
  // Sum the occupied cells 
  int newsize = mark.sum();
  
  NV.resize(newsize,V.cols());
  I.resize(V.rows(),1);
  
  // Do a pass on the marked vector and remove the unreferenced vertices
  int count = 0;
  for(int i=0;i<mark.rows();++i)
  {
    if (mark(i) == 1)
    {
      NV.row(count) = V.row(i);
      I(i) = count;
      count++;
    }
    else
    {
      I(i) = -1;
    }
  }

  NF.resize(F.rows(),F.cols());

  // Apply I on F
  for (int i=0; i<F.rows(); ++i)
  {
    Eigen::RowVectorXi t(F.cols());
    for (int j=0; j<F.cols(); ++j)
      t(j) = I(F(i,j));

    NF.row(i) = t;
  }
}

// template <typename T, typename S>
//  IGL_INLINE void removeUnreferenced(
//    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
//    const vector<vector<S> > &F,
//    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV,
//    const vector<vector<S> > &NF,
//    Eigen::Matrix<S, Eigen::Dynamic, 1> &I)
//{
//  // Mark referenced vertices
//  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> mark = Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>::Zero(V.rows(),1);
//  
//  for(int i=0; i<F.size(); ++i)
//  {
//    for(int j=0; j<F[i].size(); ++j)
//    {
//      if (F[i][j] != -1)
//        mark(F[i][j],0) = 1;
//    }
//  }
//  
//  // Sum the occupied cells 
//  int newsize = mark.sum();
//  
//  NV = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(newsize,V.cols());
//  NF.clear();
//  I  = Eigen::Matrix<S, Eigen::Dynamic, 1>(V.rows(),1);
//  
//  // Do a pass on the marked vector and remove the unreferenced vertices
//  int count = 0;
//  for(int i=0;i<mark.rows();++i)
//  {
//    if (mark(i) == 1)
//    {
//      NV.row(count) = V.row(i);
//      I(i) = count;
//      count++;
//    }
//    else
//    {
//      I(i) = -1;
//    }
//  }
//  
//  // Apply I on F
//  for (int i=0; i<F.size(); ++i)
//  {
//    vector<S> t;
//    for (int j=0; j<F[i].size(); ++j)
//      t.push_back(I[F[i][j]]);
//
//    vector<S> t_copy = t;
//    typename std::vector<S>::iterator it;
//    it = std::unique (t_copy.begin(), t_copy.end()); 
//    t_copy.resize( std::distance(t_copy.begin(),it) );
//
//    if (t_copy.size() > 2)
//      NF.push_back(t);
//  }
//}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
