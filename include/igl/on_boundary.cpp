// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "on_boundary.h"

// IGL includes
#include "sort.h"
#include "face_occurences.h"

// STL includes

template <typename IntegerT>
IGL_INLINE void igl::on_boundary(
  const std::vector<std::vector<IntegerT> > & T,
  std::vector<bool> & I,
  std::vector<std::vector<bool> > & C)
{
  using namespace std;
  using namespace igl;

  // Get a list of all faces
  vector<vector<IntegerT> > F(T.size()*4,vector<IntegerT>(3));
  // Gather faces, loop over tets
  for(int i = 0; i< (int)T.size();i++)
  {
    assert(T[i].size() == 4);
    // get face in correct order
    F[i*4+0][0] = T[i][1];
    F[i*4+0][1] = T[i][3];
    F[i*4+0][2] = T[i][2];
    // get face in correct order
    F[i*4+1][0] = T[i][0];
    F[i*4+1][1] = T[i][2];
    F[i*4+1][2] = T[i][3];
    // get face in correct order
    F[i*4+2][0] = T[i][0];
    F[i*4+2][1] = T[i][3];
    F[i*4+2][2] = T[i][1];
    // get face in correct order
    F[i*4+3][0] = T[i][0];
    F[i*4+3][1] = T[i][1];
    F[i*4+3][2] = T[i][2];
  }
  // Counts
  vector<int> FC;
  face_occurences(F,FC);
  C.resize(T.size(),vector<bool>(4));
  I.resize(T.size(),false);
  for(int i = 0; i< (int)T.size();i++)
  {
    for(int j = 0;j<4;j++)
    {
      assert(FC[i*4+j] == 2 || FC[i*4+j] == 1);
      C[i][j] = FC[i*4+j]==1;
      // if any are on boundary set to true
      I[i] = I[i] || C[i][j];
    }
  }


}

#ifndef IGL_NO_EIGEN
#include "list_to_matrix.h"
#include "matrix_to_list.h"

template <typename DerivedT, typename DerivedI, typename DerivedC>
IGL_INLINE void igl::on_boundary(
  const Eigen::PlainObjectBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedI>& I,
  Eigen::PlainObjectBase<DerivedC>& C)
{
  assert(T.cols() == 0 || T.cols() == 4);
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  // Cop out: use vector of vectors version
  vector<vector<typename Eigen::PlainObjectBase<DerivedT>::Scalar> > vT;
  matrix_to_list(T,vT);
  vector<bool> vI;
  vector<vector<bool> > vC;
  on_boundary(vT,vI,vC);
  list_to_matrix(vI,I);
  list_to_matrix(vC,C);
}
#endif


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::on_boundary<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif


