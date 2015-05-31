// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "upsample.h"

#include "triangle_triangle_adjacency.h"


template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedNV,
  typename DerivedNF>
IGL_INLINE void igl::upsample(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedNV>& NV,
  Eigen::PlainObjectBase<DerivedNF>& NF)
{
  // Use "in place" wrapper instead
  assert(&V != &NV);
  assert(&F != &NF);
  using namespace std;
  using namespace Eigen;

  Eigen::Matrix<
    typename DerivedF::Scalar,Eigen::Dynamic,Eigen::Dynamic>
    FF,FFi;
  triangle_triangle_adjacency(F,FF,FFi);

  // TODO: Cache optimization missing from here, it is a mess

  // Compute the number and positions of the vertices to insert (on edges)
  Eigen::MatrixXi NI = Eigen::MatrixXi::Constant(FF.rows(),FF.cols(),-1);
  int counter = 0;

  for(int i=0;i<FF.rows();++i)
  {
    for(int j=0;j<3;++j)
    {
      if(NI(i,j) == -1)
      {
        NI(i,j) = counter;
        if (FF(i,j) != -1) // If it is not a border
          NI(FF(i,j),FFi(i,j)) = counter;
        ++counter;
      }
    }
  }

  int n_odd = V.rows();
  int n_even = counter;

  // Preallocate NV and NF
  NV.resize(V.rows()+n_even,V.cols());
  NF.resize(F.rows()*4,3);

  // Fill the odd vertices position
  NV.block(0,0,V.rows(),V.cols()) = V;

  // Fill the even vertices position
  for(int i=0;i<FF.rows();++i)
  {
    for(int j=0;j<3;++j)
    {
      NV.row(NI(i,j) + n_odd) = 0.5 * V.row(F(i,j)) + 0.5 * V.row(F(i,(j+1)%3));
    }
  }

  // Build the new topology (Every face is replaced by four)
  for(int i=0; i<F.rows();++i)
  {
    VectorXi VI(6);
    VI << F(i,0), F(i,1), F(i,2), NI(i,0) + n_odd, NI(i,1) + n_odd, NI(i,2) + n_odd;

    VectorXi f0(3), f1(3), f2(3), f3(3);
    f0 << VI(0), VI(3), VI(5);
    f1 << VI(1), VI(4), VI(3);
    f2 << VI(3), VI(4), VI(5);
    f3 << VI(4), VI(2), VI(5);

    NF.row((i*4)+0) = f0;
    NF.row((i*4)+1) = f1;
    NF.row((i*4)+2) = f2;
    NF.row((i*4)+3) = f3;
  }

}

template <
  typename MatV,
  typename MatF>
IGL_INLINE void igl::upsample(
  MatV& V,
  MatF& F)
{
  const MatV V_copy = V;
  const MatF F_copy = F;
  return upsample(V_copy,F_copy,V,F);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::upsample<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::Matrix<double, -1, -1, 0, -1, -1>&, Eigen::Matrix<int, -1, -1, 0, -1, -1>&);
#endif
