// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unique_simplices.h"
#include "sort.h"
#include "unique.h"
#include "get_seconds.h"

template <
  typename DerivedF,
  typename DerivedFF,
  typename DerivedIA,
  typename DerivedIC>
IGL_INLINE void igl::unique_simplices(
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedFF>& FF,
  Eigen::PlainObjectBase<DerivedIA>& IA,
  Eigen::PlainObjectBase<DerivedIC>& IC)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  // Sort each face
  MatrixXi sortF, unusedI;
  igl::sort(F,2,true,sortF,unusedI);
  // Find unique faces
  MatrixXi C;
  igl::unique_rows(sortF,C,IA,IC);
  FF.resize(IA.size(),F.cols());
  const size_t mff = FF.rows();
  // Minimum number of iterms per openmp thread
  #ifndef IGL_OMP_MIN_VALUE
  #  define IGL_OMP_MIN_VALUE 1000
  #endif
  #pragma omp parallel for if (mff>IGL_OMP_MIN_VALUE)
  // Copy into output
  for(size_t i = 0;i<mff;i++)
  {
    FF.row(i) = F.row(IA(i));
  }
}

template <
  typename DerivedF,
  typename DerivedFF>
IGL_INLINE void igl::unique_simplices(
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedFF>& FF)
{
  Eigen::VectorXi IA,IC;
  return unique_simplices(F,FF,IA,IC);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instanciations
template void igl::unique_simplices<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::unique_simplices<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template void igl::unique_simplices<Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
#endif
