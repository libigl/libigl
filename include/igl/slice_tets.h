// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SLICE_TETS_H
#define IGL_SLICE_TETS_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace igl
{
  // SLICE_TETS Slice through a tet mesh (V,T) along a given plane (via its
  // implicit equation).
  //
  // Inputs:
  //   V  #V by 3 list of tet mesh vertices
  //   T  #T by 4 list of tet indices into V 
  //   plane  list of 4 coefficients in the plane equation: [x y z 1]'*plane = 0
  //   Optional:
  //     'Manifold' followed by whether to stitch together triangles into a
  //       manifold mesh {true}: results in more compact U but slightly slower.
  // Outputs:
  //   U  #U by 3 list of triangle mesh vertices along slice
  //   G  #G by 3 list of triangles indices into U
  //   J  #G list of indices into T revealing from which tet each faces comes
  //   BC  #U by #V list of barycentric coordinates (or more generally: linear
  //     interpolation coordinates) so that U = BC*V
  // 
  template <
    typename DerivedV, 
    typename DerivedT, 
    typename Derivedplane,
    typename DerivedU,
    typename DerivedG,
    typename DerivedJ,
    typename BCType>
  IGL_INLINE void slice_tets(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const Eigen::PlainObjectBase<Derivedplane> & plane,
    Eigen::PlainObjectBase<DerivedU>& U,
    Eigen::PlainObjectBase<DerivedG>& G,
    Eigen::PlainObjectBase<DerivedJ>& J,
    Eigen::SparseMatrix<BCType> & BC);
}

#ifndef IGL_STATIC_LIBRARY
#  include "slice_tets.cpp"
#endif

#endif


