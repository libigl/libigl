// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef CROUZEIX_RAVIART_MASSMATRIX_H
#define CROUZEIX_RAVIART_MASSMATRIX_H
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // CROUZEIX_RAVIART_MASSMATRIX Compute the Crouzeix-Raviart mass matrix where
  // M(e,e) is just the sum of the areas of the triangles on either side of an
  // edge e.
  //
  // See for example "Discrete Quadratic Curvature Energies" [Wardetzky, Bergou,
  // Harmon, Zorin, Grinspun 2007]
  //
  // Templates:
  //   MT  type of eigen sparse matrix for M (e.g. double for
  //     SparseMatrix<double>)
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   DerivedE  derived type of eigen matrix for E (e.g. derived from
  //     MatrixXi)
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   M  #E by #E edge-based diagonal mass matrix
  //   E  #E by 2 list of edges
  //   EMAP  #F*3 list of indices mapping allE to E
  //
  //
  template <typename MT, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP>
  void crouzeix_raviart_massmatrix(
      const Eigen::PlainObjectBase<DerivedV> & V, 
      const Eigen::PlainObjectBase<DerivedF> & F, 
      Eigen::SparseMatrix<MT> & M,
      Eigen::PlainObjectBase<DerivedE> & E,
      Eigen::PlainObjectBase<DerivedEMAP> & EMAP);
  // wrapper if E and EMAP are already computed (better match!)
  template <typename MT, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP>
  void crouzeix_raviart_massmatrix(
      const Eigen::PlainObjectBase<DerivedV> & V, 
      const Eigen::PlainObjectBase<DerivedF> & F, 
      const Eigen::PlainObjectBase<DerivedE> & E,
      const Eigen::PlainObjectBase<DerivedEMAP> & EMAP,
      Eigen::SparseMatrix<MT> & M);
}
#ifndef IGL_STATIC_LIBRARY
#  include "crouzeix_raviart_massmatrix.cpp"
#endif
  
#endif
