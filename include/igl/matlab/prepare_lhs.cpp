// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "prepare_lhs.h"
#include <algorithm>
template <typename DerivedV>
IGL_INLINE void igl::matlab::prepare_lhs_double(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  using namespace std;
  using namespace Eigen;
  plhs[0] = mxCreateDoubleMatrix(V.rows(),V.cols(), mxREAL);
  double * Vp = mxGetPr(plhs[0]);

  typedef typename DerivedV::Scalar Scalar;
  const Scalar * V_data;
  Matrix<Scalar, DerivedV::ColsAtCompileTime, DerivedV::RowsAtCompileTime, RowMajor> VT;
  if(DerivedV::IsRowMajor)
  {
    VT = V.transpose();
    V_data = VT.data();
  }else
  {
    V_data = V.data();
  }

  copy(V_data,V_data+V.size(),Vp);
}

template <typename DerivedV>
IGL_INLINE void igl::matlab::prepare_lhs_logical(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  using namespace std;
  using namespace Eigen;
  plhs[0] = mxCreateLogicalMatrix(V.rows(),V.cols());
  mxLogical * Vp = static_cast<mxLogical*>(mxGetData(plhs[0]));

  typedef typename DerivedV::Scalar Scalar;
  const Scalar * V_data;
  Matrix<Scalar, DerivedV::ColsAtCompileTime, DerivedV::RowsAtCompileTime, RowMajor> VT;
  if(DerivedV::IsRowMajor)
  {
    VT = V.transpose();
    V_data = VT.data();
  }else
  {
    V_data = V.data();
  }

  copy(V_data,V_data+V.size(),Vp);
}

template <typename DerivedV>
IGL_INLINE void igl::matlab::prepare_lhs_index(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  // Treat indices as reals
  
  const auto Vd = (V.template cast<double>().array()+1).eval();
  return prepare_lhs_double(Vd,plhs);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::matlab::prepare_lhs_index<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
template void igl::matlab::prepare_lhs_index<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
template void igl::matlab::prepare_lhs_double<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, mxArray_tag**);
template void igl::matlab::prepare_lhs_index<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, mxArray_tag**);
template void igl::matlab::prepare_lhs_logical<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
template void igl::matlab::prepare_lhs_double<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
#endif
