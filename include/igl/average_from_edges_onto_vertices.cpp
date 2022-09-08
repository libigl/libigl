// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "average_from_edges_onto_vertices.h"

template<typename DerivedF,typename DerivedE,typename DerivedoE,
typename DeriveduE,typename DeriveduV>
IGL_INLINE void
igl::average_from_edges_onto_vertices(
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedE> &E,
  const Eigen::MatrixBase<DerivedoE> &oE,
  const Eigen::MatrixBase<DeriveduE> &uE,
  Eigen::PlainObjectBase<DeriveduV> &uV)
{
  using Scalar = typename DeriveduE::Scalar;
  using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using Int = typename DerivedF::Scalar;
  
  assert(E.rows()==F.rows() && "E does not match dimensions of F.");
  assert(oE.rows()==F.rows() && "oE does not match dimensions of F.");
  assert(E.cols()==3 && F.cols()==3 && oE.cols()==3 &&
   "This method is for triangle meshes.");
  
  const Int n = F.maxCoeff()+1;
  
  VecX edgesPerVertex(n);
  edgesPerVertex.setZero();
  uV.resize(n,1);
  uV.setZero();
  
  for(Eigen::Index i=0; i<F.rows(); ++i) {
    for(int j=0; j<3; ++j) {
      if(oE(i,j)<0) {
        continue;
      }
      const Int e = E(i,j);
      const Int vi=F(i,(j+1)%3), vj=F(i,(j+2)%3);
      
      //Count vertex valence
      ++edgesPerVertex(vi);
      ++edgesPerVertex(vj);
      
      //Average uE value onto vertices
      uV(vi) += uE(e);
      uV(vj) += uE(e);
    }
  }
  
  //Divide by valence
  for(Int i=0; i<n; ++i) {
    const Scalar valence = edgesPerVertex(i);
    if(valence>0) {
      uV(i) /= valence;
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::average_from_edges_onto_vertices<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::PartialReduxExpr<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::internal::member_norm<double>, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::PartialReduxExpr<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::internal::member_norm<double>, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::average_from_edges_onto_vertices<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::average_from_edges_onto_vertices<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
