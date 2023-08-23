// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich <michaelrabinovich27@gmail.com@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "euler_characteristic.h"

#include "unique_edge_map.h"
#include <cassert>

template <typename DerivedF>
IGL_INLINE int igl::euler_characteristic(
  const Eigen::MatrixBase<DerivedF> & F)
{
  const int nf = F.rows();
  const int nv = F.maxCoeff()+1;
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,2> E;
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> EMAP;
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,2> uE;
  unique_edge_map(F,E,uE,EMAP);
  // The input is assumed manifold so boundary edges are those with exactly one
  // incident face
  Eigen::VectorXi count = Eigen::VectorXi::Zero(uE.rows(),1);
  for(int e = 0;e<EMAP.size();e++)
  {
    count(EMAP(e))++;
  }

  std::vector<Eigen::Triplet<int> > IJV;
  for(int u = 0;u<uE.rows();u++)
  {
    if(count(u) == 1)
    {
      IJV.emplace_back(uE(u,0),uE(u,1),1);
      IJV.emplace_back(uE(u,1),uE(u,0),1);
    }
    assert(count(u) == 2);
  }
  Eigen::SparseMatrix<int> A(nv,nv);
  A.setFromTriplets(IJV.begin(),IJV.end());
  Eigen::VectorXi _C,_K;
  // Number of boundary connected components
  const int nb = connected_components(A,_C,_K);
  std::cout<<"nb: "<<nb<<std::endl;

  const int ne = E.rows();
  return nv - ne + nf + nb;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::euler_characteristic<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
#endif
