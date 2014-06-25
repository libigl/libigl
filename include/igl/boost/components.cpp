// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "components.h"
#include <igl/adjacency_matrix.h>

//#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <iostream>
#include <vector>
#include <cassert>

template <typename AScalar, typename DerivedC>
IGL_INLINE void igl::components(
  const Eigen::SparseMatrix<AScalar> & A,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  assert(A.rows() == A.cols());
  using namespace Eigen;
  // THIS IS DENSE:
  //boost::adjacency_matrix<boost::undirectedS> bA(A.rows());
  boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> bA(A.rows());
  for(int j=0; j<A.outerSize();j++)
  {
    // Iterate over inside
    for(typename SparseMatrix<AScalar>::InnerIterator it (A,j); it; ++it)
    {
      if(0 != it.value())
      {
        boost::add_edge(it.row(),it.col(),bA);
      }
    }
  }
  C.resize(A.rows(),1);
  boost::connected_components(bA,C.data());
}

template <typename DerivedF, typename DerivedC>
IGL_INLINE void igl::components(
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  Eigen::SparseMatrix<typename DerivedC::Scalar> A;
  igl::adjacency_matrix(F,A);
  return components(A,C);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::components<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
