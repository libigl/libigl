// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "average_onto_faces.h"

template <typename DerivedF, typename DerivedS, typename DerivedSF>
IGL_INLINE void average_onto_faces(
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedS> & S,
  Eigen::PlainObjectBase<DerivedSF> & SF)
{
  SF.setConstant(F.rows(),S.cols(),0);
  for (int i = 0; i <F.rows(); ++i)
    for (int j = 0; j<F.cols(); ++j)
      SF.row(i) += S.row(F(i,j));
  SF.array() /= F.cols();
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
