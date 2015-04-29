// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "crouzeix_raviart_massmatrix.h"
#include "unique_simplices.h"
#include "all_edges.h"

#include "is_edge_manifold.h"
#include "doublearea.h"

#include <cassert>
#include <vector>

template <typename MT, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP>
void igl::crouzeix_raviart_massmatrix(
    const Eigen::PlainObjectBase<DerivedV> & V, 
    const Eigen::PlainObjectBase<DerivedF> & F, 
    Eigen::SparseMatrix<MT> & M,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DerivedEMAP> & EMAP)
{
  // All occurances of directed edges
  Eigen::MatrixXi allE;
  all_edges(F,allE);
  Eigen::VectorXi _1;
  unique_simplices(allE,E,_1,EMAP);
  return crouzeix_raviart_massmatrix(V,F,E,EMAP,M);
}

template <typename MT, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP>
void igl::crouzeix_raviart_massmatrix(
    const Eigen::PlainObjectBase<DerivedV> & V, 
    const Eigen::PlainObjectBase<DerivedF> & F, 
    const Eigen::PlainObjectBase<DerivedE> & E,
    const Eigen::PlainObjectBase<DerivedEMAP> & EMAP,
    Eigen::SparseMatrix<MT> & M)
{
  using namespace Eigen;
  using namespace std;
  assert(F.cols() == 3);
  // Mesh should be edge-manifold
  assert(is_edge_manifold(V,F));
  // number of elements (triangles)
  int m = F.rows();
  // Get triangle areas
  VectorXd TA;
  doublearea(V,F,TA);
  vector<Triplet<MT> > MIJV(3*m);
  for(int f = 0;f<m;f++)
  {
    for(int c = 0;c<3;c++)
    {
      MIJV[f+m*c] = Triplet<MT>(EMAP(f+m*c),EMAP(f+m*c),TA(f)*0.5);
    }
  }
  M.resize(E.rows(),E.rows());
  M.setFromTriplets(MIJV.begin(),MIJV.end());
}

#ifdef IGL_STATIC_LIBRARY
// Explicit instanciation
template void igl::crouzeix_raviart_massmatrix<double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int>&);
template void igl::crouzeix_raviart_massmatrix<float, Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<float, 0, int>&);
#endif
