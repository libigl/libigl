// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "split_nonmanifold.h"
#include "connected_components.h"
#include "remove_unreferenced.h"
#include "find.h"
#include "ismember_rows.h"

template <
  typename DerivedF,
  typename DerivedSF,
  typename DerivedSVI
  >
IGL_INLINE void igl::split_nonmanifold(
  const Eigen::MatrixBase<DerivedF> & F,
  Eigen::PlainObjectBase <DerivedSF> & SF,
  Eigen::PlainObjectBase <DerivedSVI> & SVI)
{
  // Number of faces
  const int m = F.rows();
  // For moment aact like everything will be split
  SF.resize(m,3);
  {
    int k =0;
    for(int j = 0;j<3;j++)
    {
      for(int i = 0;i<m;i++)
      {
        SF(i,j) = k++;
      }
    }
  }
  // Edges in SF
  Eigen::MatrixXi E(m*3,2);
  for(int i = 0;i<m;i++)
  {
    E.row(i+0*m) << SF(i,1),SF(i,2);
    E.row(i+1*m) << SF(i,2),SF(i,0);
    E.row(i+2*m) << SF(i,0),SF(i,1);
  }
  // Reindex E by F
  Eigen::MatrixXi FE(E.rows(),E.cols());
  for(int i = 0;i<E.rows();i++)
  {
    for(int j = 0;j<2;j++)
    {
      const int fi = E(i,j) % m;
      const int fj = E(i,j) / m;
      FE(i,j) = F(fi,fj);
    }
  }
  // Flip orientation
  Eigen::MatrixXi FE_flip = FE.rowwise().reverse();
  // Find which exist in both directions


  Eigen::Array<bool,Eigen::Dynamic,1> I;
  Eigen::VectorXi J;
  igl::ismember_rows(FE,FE_flip,I,J);
  // Just keep those find
  const auto II = igl::find(I);
  Eigen::MatrixXi EI = E(II,Eigen::all);
  Eigen::VectorXi JI = J(II);
  Eigen::MatrixXi EJI = E(JI,Eigen::all);

  Eigen::MatrixXi EJI_flip = EJI.rowwise().reverse();
  // Build adjacency matrix
  std::vector<Eigen::Triplet<bool> > Aijv; 
  Aijv.reserve(EI.size());
  for(int i = 0;i<EI.rows();i++)
  {
    for(int j = 0;j<2;j++)
    {
      Aijv.emplace_back( 
        EI(i,j), 
        EJI_flip(i,j), 
        true);
    }
  }
  // Build A to contain off-diagonals only if both directions are present
  Eigen::SparseMatrix<bool> A1(m*3,m*3);
  A1.setFromTriplets(Aijv.begin(),Aijv.end());
  // For some reason I can't write `A = A1 && A1.transpose();`
  Eigen::SparseMatrix<bool> A1T = A1.transpose();
  Eigen::SparseMatrix<bool> A = A1 && A1T;


  Eigen::VectorXi K;
  {
    Eigen::VectorXi _;
    igl::connected_components(A,K,_);
  }


  // Remap by components
  for(int j = 0;j<3;j++)
  {
    for(int i = 0;i<m;i++)
    {
      SF(i,j) = K(SF(i,j));
    }
  }
  
  // Initial mapping
  Eigen::VectorXi SVI0(m*3);
  {
    int k =0;
    for(int j = 0;j<3;j++)
    {
      for(int i = 0;i<m;i++)
      {
        SVI0(k++) = F(i,j);
      }
    }
  }
  assert(K.size() == m*3);
  // Scatter via K
  // SVI1(K) = SVI(K);
  Eigen::VectorXi SVI1(m*3);
  SVI1(K) = SVI0;

  {
    Eigen::VectorXi _,J;
    igl::remove_unreferenced(SF.maxCoeff()+1,SF,_,J);
    // Remap by J
    for(int j = 0;j<3;j++)
    {
      for(int i = 0;i<m;i++)
      {
        SF(i,j) = J(SF(i,j));
      }
    }
    SVI = SVI1(J);
  }

}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedSV,
  typename DerivedSF,
  typename DerivedSVI
  >
IGL_INLINE void igl::split_nonmanifold(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  Eigen::PlainObjectBase <DerivedSV> & SV,
  Eigen::PlainObjectBase <DerivedSF> & SF,
  Eigen::PlainObjectBase <DerivedSVI> & SVI)
{
  igl::split_nonmanifold(F,SF,SVI);
  SV = V(SVI.derived(),Eigen::all);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::split_nonmanifold<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
