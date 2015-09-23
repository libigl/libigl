// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "remove_unreferenced.h"
#include "slice.h"
#include <algorithm>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedNV,
  typename DerivedNF,
  typename DerivedI>
IGL_INLINE void igl::remove_unreferenced(
  const Eigen::PlainObjectBase<DerivedV> &V,
  const Eigen::PlainObjectBase<DerivedF> &F,
  Eigen::PlainObjectBase<DerivedNV> &NV,
  Eigen::PlainObjectBase<DerivedNF> &NF,
  Eigen::PlainObjectBase<DerivedI> &I)
{
  Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,1> J;
  remove_unreferenced(V,F,NV,NF,I,J);
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedNV,
  typename DerivedNF,
  typename DerivedI,
  typename DerivedJ>
IGL_INLINE void igl::remove_unreferenced(
  const Eigen::PlainObjectBase<DerivedV> &V,
  const Eigen::PlainObjectBase<DerivedF> &F,
  Eigen::PlainObjectBase<DerivedNV> &NV,
  Eigen::PlainObjectBase<DerivedNF> &NF,
  Eigen::PlainObjectBase<DerivedI> &I,
  Eigen::PlainObjectBase<DerivedJ> &J)
{
  using namespace std;
  const size_t n = V.rows();
  remove_unreferenced(n,F,I,J);
  NF = F;
  for_each(NF.data(),NF.data()+NF.size(),[&I](int & a){a=I(a);});
  slice(V,J,1,NV);
}

template <
  typename DerivedF,
  typename DerivedI,
  typename DerivedJ>
IGL_INLINE void igl::remove_unreferenced(
  const size_t n,
  const Eigen::PlainObjectBase<DerivedF> &F,
  Eigen::PlainObjectBase<DerivedI> &I,
  Eigen::PlainObjectBase<DerivedJ> &J)
{
  // Mark referenced vertices
  typedef Eigen::Matrix<bool,Eigen::Dynamic,1> MatrixXb;
  MatrixXb mark = MatrixXb::Zero(n,1);
  for(int i=0; i<F.rows(); ++i)
  {
    for(int j=0; j<F.cols(); ++j)
    {
      if (F(i,j) != -1)
      {
        mark(F(i,j)) = 1;
      }
    }
  }

  // Sum the occupied cells
  int newsize = mark.count();

  I.resize(n,1);
  J.resize(newsize,1);

  // Do a pass on the marked vector and remove the unreferenced vertices
  int count = 0;
  for(int i=0;i<mark.rows();++i)
  {
    if (mark(i) == 1)
    {
      I(i) = count;
      J(count) = i;
      count++;
    }
    else
    {
      I(i) = -1;
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::remove_unreferenced<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::remove_unreferenced<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::remove_unreferenced<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
