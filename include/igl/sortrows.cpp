// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sortrows.h"
#include "get_seconds.h"

#include "SortableRow.h"
#include "sort.h"
#include "colon.h"
#include "IndexComparison.h"

#include <vector>

// Obsolete slower version converst to vector
//template <typename DerivedX, typename DerivedIX>
//IGL_INLINE void igl::sortrows(
//  const Eigen::PlainObjectBase<DerivedX>& X,
//  const bool ascending,
//  Eigen::PlainObjectBase<DerivedX>& Y,
//  Eigen::PlainObjectBase<DerivedIX>& IX)
//{
//  using namespace std;
//  using namespace Eigen;
//  typedef Eigen::Matrix<typename DerivedX::Scalar, Eigen::Dynamic, 1> RowVector;
//  vector<SortableRow<RowVector> > rows;
//  rows.resize(X.rows());
//  // Loop over rows
//  for(int i = 0;i<X.rows();i++)
//  {
//    RowVector ri = X.row(i);
//    rows[i] = SortableRow<RowVector>(ri);
//  }
//  vector<SortableRow<RowVector> > sorted;
//  std::vector<size_t> index_map;
//  // Perform sort on rows
//  igl::sort(rows,ascending,sorted,index_map);
//  // Resize output
//  Y.resize(X.rows(),X.cols());
//  IX.resize(X.rows(),1);
//  // Convert to eigen
//  for(int i = 0;i<X.rows();i++)
//  {
//    Y.row(i) = sorted[i].data;
//    IX(i,0) = index_map[i];
//  }
//}

template <typename DerivedX, typename DerivedIX>
IGL_INLINE void igl::sortrows(
  const Eigen::PlainObjectBase<DerivedX>& X,
  const bool ascending,
  Eigen::PlainObjectBase<DerivedX>& Y,
  Eigen::PlainObjectBase<DerivedIX>& IX)
{
  // This is already 2x faster than matlab's builtin `sortrows`. I have tried
  // implementing a "multiple-pass" sort on each column, but see no performance
  // improvement.
  using namespace std;
  using namespace Eigen;
  // Resize output
  Y.resize(X.rows(),X.cols());
  IX.resize(X.rows(),1);
  for(int i = 0;i<X.rows();i++)
  {
    IX(i) = i;
  }
  std::sort(
    IX.data(),
    IX.data()+IX.size(),
    igl::IndexRowLessThan<const Eigen::PlainObjectBase<DerivedX> & >(X));
  // if not ascending then reverse
  if(!ascending)
  {
    std::reverse(IX.data(),IX.data()+IX.size());
  }
  for(int i = 0;i<X.rows();i++)
  {
    Y.row(i) = X.row(IX(i));
  }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::sortrows<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sortrows<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sortrows<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::sortrows<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sortrows<Eigen::Matrix<double, -1, 2, 0, -1, 2>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sortrows<Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
