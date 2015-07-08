// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "cdt.h"
#include "../bounding_box.h"
#include "../writeOBJ.h"

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedTV, 
  typename DerivedTT, 
  typename DerivedTF>
IGL_INLINE bool igl::tetgen::cdt(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const CDTParam & param,
  Eigen::PlainObjectBase<DerivedTV>& TV,
  Eigen::PlainObjectBase<DerivedTT>& TT,
  Eigen::PlainObjectBase<DerivedTF>& TF)
{
  using namespace Eigen;
  using namespace std;
  typedef Eigen::PlainObjectBase<DerivedV> MatrixXS;
  typedef Eigen::PlainObjectBase<DerivedF> MatrixXI;
  // Effective input mesh
  MatrixXS U;
  MatrixXI G;
  if(param.use_bounding_box)
  {
    // Construct bounding box mesh
    MatrixXS BV;
    MatrixXI BF;
    bounding_box(V,BV,BF);
    // scale bounding box
    const RowVector3d mid = 
     (BV.colwise().minCoeff() + BV.colwise().maxCoeff()).eval()*0.5;
    BV.rowwise() -= mid;
    assert(param.bounding_box_scale >= 1.);
    BV.array() *= param.bounding_box_scale;
    BV.rowwise() += mid;
    // Append bounding box to mesh
    U.resize(V.rows()+BV.rows(),V.cols());
    U<<V,BV;
    BF.array() += V.rows();
    G.resize(F.rows()+BF.rows(),F.cols());
    G<<F,BF;
  }else
  {
    // needless copies
    U = V;
    G = F;
  }
  // effective flags;
  string flags = param.flags + (param.use_bounding_box ? "" : "c");
  writeOBJ("UG.obj",U,G);
  return tetrahedralize(U,G,flags,TV,TT,TF);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
