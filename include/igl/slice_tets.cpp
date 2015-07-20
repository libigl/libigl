// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "slice_tets.h"
#include <igl/sort.h>
#include <igl/cat.h>
#include <igl/per_face_normals.h>
#include <cassert>
#include <algorithm>
#include <vector>

template <
  typename DerivedV, 
  typename DerivedT, 
  typename Derivedplane,
  typename DerivedU,
  typename DerivedG,
  typename DerivedJ,
  typename BCType>
IGL_INLINE void igl::slice_tets(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T,
  const Eigen::PlainObjectBase<Derivedplane> & plane,
  Eigen::PlainObjectBase<DerivedU>& U,
  Eigen::PlainObjectBase<DerivedG>& G,
  Eigen::PlainObjectBase<DerivedJ>& J,
  Eigen::SparseMatrix<BCType> & BC)
{
  using namespace Eigen;
  using namespace std;
  assert(V.cols() == 3 && "V should be #V by 3");
  assert(T.cols() == 4 && "T should be #T by 4");
  assert(plane.size() == 4 && "Plane equation should be 4 coefficients");

  // number of tets
  const size_t m = T.rows();

  typedef typename DerivedV::Scalar Scalar;
  typedef typename DerivedT::Scalar Index;
  typedef Matrix<Scalar,Dynamic,1> VectorXS;
  typedef Matrix<Scalar,Dynamic,4> MatrixX4S;
  typedef Matrix<Scalar,Dynamic,3> MatrixX3S;
  typedef Matrix<Scalar,Dynamic,2> MatrixX2S;
  typedef Matrix<Index,Dynamic,4> MatrixX4I;
  typedef Matrix<Index,Dynamic,3> MatrixX3I;
  typedef Matrix<Index,Dynamic,1> VectorXI;
  typedef Matrix<bool,Dynamic,1> VectorXb;
  
  // Value of plane's implicit function at all vertices
  VectorXS IV = 
    (V.col(0)*plane(0) + 
     V.col(1)*plane(1) + 
     V.col(2)*plane(2)).array()
    + plane(3);
  MatrixX4S IT(m,4);
  for(size_t t = 0;t<m;t++)
  {
    for(size_t c = 0;c<4;c++)
    {
      IT(t,c) = IV(T(t,c));
    }
  }

  const auto & extract_rows = [](
    const PlainObjectBase<DerivedT> & T,
    const MatrixX4S & IT,
    const VectorXb & I,
    MatrixX4I  & TI,
    MatrixX4S & ITI,
    VectorXI & JI)
  {
    const Index num_I = std::count(I.data(),I.data()+I.size(),true);
    TI.resize(num_I,4);
    ITI.resize(num_I,4);
    JI.resize(num_I,1);
    {
      size_t k = 0;
      for(size_t t = 0;t<(size_t)T.rows();t++)
      {
        if(I(t))
        {
          TI.row(k) = T.row(t);
          ITI.row(k) = IT.row(t);
          JI(k) = t;
          k++;
        }
      }
      assert(k == num_I);
    }
  };

  VectorXb I13 = (IT.array()<0).rowwise().count()==1;
  VectorXb I31 = (IT.array()>0).rowwise().count()==1;
  VectorXb I22 = (IT.array()<0).rowwise().count()==2;
  MatrixX4I T13,T31,T22;
  MatrixX4S IT13,IT31,IT22;
  VectorXI J13,J31,J22;
  extract_rows(T,IT,I13,T13,IT13,J13);
  extract_rows(T,IT,I31,T31,IT31,J31);
  extract_rows(T,IT,I22,T22,IT22,J22);

  const auto & apply_sort = [] (
     const MatrixX4I & T, 
     const MatrixX4I & sJ, 
     MatrixX4I & sT)
  {
    sT.resize(T.rows(),4);
    for(size_t t = 0;t<(size_t)T.rows();t++)
    {
      for(size_t c = 0;c<4;c++)
      {
        sT(t,c) = T(t,sJ(t,c));
      }
    }
  };

  const auto & one_below = [&V,&apply_sort](
    const MatrixX4I & T,
    const MatrixX4S & IT,
    MatrixX3I & G,
    SparseMatrix<BCType> & BC)
  {
    // Number of tets
    const size_t m = T.rows();
    MatrixX4S sIT;
    MatrixX4I sJ;
    sort(IT,2,true,sIT,sJ);
    MatrixX4I sT;
    apply_sort(T,sJ,sT);
    MatrixX3S lambda = 
      sIT.rightCols(3).array() /
      (sIT.rightCols(3).colwise()-sIT.col(0)).array();
    vector<Triplet<BCType> > IJV;
    IJV.reserve(m*3*2);
    for(size_t c = 0;c<3;c++)
    {
      for(size_t t = 0;t<(size_t)m;t++)
      {
        IJV.push_back(Triplet<BCType>(c*m+t,  sT(t,0),  lambda(t,c)));
        IJV.push_back(Triplet<BCType>(c*m+t,sT(t,c+1),1-lambda(t,c)));
      }
    }
    BC.resize(m*3,V.rows());
    BC.reserve(m*3*2);
    BC.setFromTriplets(IJV.begin(),IJV.end());
    G.resize(m,3);
    for(size_t c = 0;c<3;c++)
    {
      G.col(c).setLinSpaced(m,0+c*m,(m-1)+c*m);
    }
  };

  const auto & two_below = [&V,&apply_sort](
    const MatrixX4I & T,
    const MatrixX4S & IT,
    MatrixX3I & G,
    SparseMatrix<BCType> & BC)
  {
    // Number of tets
    const size_t m = T.rows();
    MatrixX4S sIT;
    MatrixX4I sJ;
    sort(IT,2,true,sIT,sJ);
    MatrixX4I sT;
    apply_sort(T,sJ,sT);
    MatrixX2S lambda = 
      sIT.rightCols(2).array() /
      (sIT.rightCols(2).colwise()-sIT.col(0)).array();
    MatrixX2S gamma = 
      sIT.rightCols(2).array() /
      (sIT.rightCols(2).colwise()-sIT.col(1)).array();
    vector<Triplet<BCType> > IJV;
    IJV.reserve(m*4*2);
    for(size_t c = 0;c<2;c++)
    {
      for(size_t t = 0;t<(size_t)m;t++)
      {
        IJV.push_back(Triplet<BCType>(0*2*m+c*m+t,  sT(t,0),  lambda(t,c)));
        IJV.push_back(Triplet<BCType>(0*2*m+c*m+t,sT(t,c+2),1-lambda(t,c)));
        IJV.push_back(Triplet<BCType>(1*2*m+c*m+t,  sT(t,1),   gamma(t,c)));
        IJV.push_back(Triplet<BCType>(1*2*m+c*m+t,sT(t,c+2),1- gamma(t,c)));
      }
    }
    BC.resize(m*4,V.rows());
    BC.reserve(m*4*2);
    BC.setFromTriplets(IJV.begin(),IJV.end());
    G.resize(2*m,3);
    G.block(0,0,m,1) = VectorXI::LinSpaced(m,0+0*m,(m-1)+0*m);
    G.block(0,1,m,1) = VectorXI::LinSpaced(m,0+1*m,(m-1)+1*m);
    G.block(0,2,m,1) = VectorXI::LinSpaced(m,0+3*m,(m-1)+3*m);
    G.block(m,0,m,1) = VectorXI::LinSpaced(m,0+0*m,(m-1)+0*m);
    G.block(m,1,m,1) = VectorXI::LinSpaced(m,0+3*m,(m-1)+3*m);
    G.block(m,2,m,1) = VectorXI::LinSpaced(m,0+2*m,(m-1)+2*m);
  };

  MatrixX3I G13,G31,G22;
  SparseMatrix<BCType> BC13,BC31,BC22;
  one_below(T13,IT13,G13,BC13);
  one_below(T31,-IT31,G31,BC31);
  two_below(T22,IT22,G22,BC22);

  BC = cat(1,cat(1,BC13,BC31),BC22);
  U = BC*V;
  G.resize(G13.rows()+G31.rows()+G22.rows(),3);
  G<<G13,(G31.array()+BC13.rows()),(G22.array()+BC13.rows()+BC31.rows());
  MatrixX3S N;
  per_face_normals(U,G,N);
  Matrix<Scalar,1,3> planeN(plane(0),plane(1),plane(2));
  VectorXb flip = (N.array().rowwise() * planeN.array()).rowwise().sum()<0;
  for(size_t g = 0;g<(size_t)G.rows();g++)
  {
    if(flip(g))
    {
      G.row(g) = G.row(g).reverse().eval();
    }
  }

  J.resize(G.rows());
  J<<J13,J31,J22,J22;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::slice_tets<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 4, 1, 0, 4, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, double>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 4, 1, 0, 4, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::SparseMatrix<double, 0, int>&);
#endif
