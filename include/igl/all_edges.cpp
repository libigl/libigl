// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "all_edges.h"

template <typename DerivedF, typename DerivedE>
IGL_INLINE void igl::all_edges(
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedE> & E)
{
  E.resize(F.rows()*F.cols(),F.cols()-1);
  switch(F.cols())
  {
    case 4:
      E.block(0*F.rows(),0,F.rows(),1) = F.col(1);
      E.block(0*F.rows(),1,F.rows(),1) = F.col(3);
      E.block(0*F.rows(),2,F.rows(),1) = F.col(2);

      E.block(1*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(1*F.rows(),1,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),2,F.rows(),1) = F.col(3);

      E.block(2*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),1,F.rows(),1) = F.col(3);
      E.block(2*F.rows(),2,F.rows(),1) = F.col(1);

      E.block(3*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(3*F.rows(),1,F.rows(),1) = F.col(1);
      E.block(3*F.rows(),2,F.rows(),1) = F.col(2);
      return;
    case 3:
      E.block(0*F.rows(),0,F.rows(),1) = F.col(1);
      E.block(0*F.rows(),1,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),0,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),1,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),1,F.rows(),1) = F.col(1);
      return;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::all_edges<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::all_edges<Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&);
template void igl::all_edges<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&);
template void igl::all_edges<Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::all_edges<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 2, 0, -1, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> >&);
#endif
