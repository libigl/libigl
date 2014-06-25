// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "vector_area_matrix.h"
#include <vector>

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

//#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>

template <typename DerivedF, typename Scalar>
IGL_INLINE void igl::vector_area_matrix(
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::SparseMatrix<Scalar>& A)
{
  using namespace Eigen;
  using namespace std;

  // number of vertices
  const int n = F.maxCoeff()+1;

	SparseMatrix<Scalar> aux (n * 2, n * 2);
	SparseMatrix<Scalar> auxT(n * 2, n * 2);

	vector<Triplet<Scalar> > auxTripletList;
	vector<Triplet<Scalar> > auxTTripletList;

  MatrixXi E;
  boundary_facets(F,E);

	for(int k = 0; k < E.rows(); k++)
  {
		int i = E(k,0);
		int j = E(k,1);
		auxTripletList.push_back(Triplet<Scalar>(i+n, j, -0.5));
		auxTripletList.push_back(Triplet<Scalar>(i, j+n, 0.5));
		auxTTripletList.push_back(Triplet<Scalar>(j, i+n, -0.5));
		auxTTripletList.push_back(Triplet<Scalar>(j+n, i, 0.5));
	}

	aux.setFromTriplets(auxTripletList.begin(), auxTripletList.end());
	auxT.setFromTriplets(auxTTripletList.begin(), auxTTripletList.end());
	A = (aux + auxT)/2;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::vector_area_matrix<Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int>&);
#endif
