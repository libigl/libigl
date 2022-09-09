// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Xiangyu Kong <xiangyu.kong@mail.utoronto.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DIRECT_DELTA_MUSH_H
#define IGL_DIRECT_DELTA_MUSH_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <vector>

namespace igl {
  // Computes Direct Delta Mesh Skinning (Variant 0) from "Direct Delta Mush
  // Skinning and Variants"
  //
  // Inputs:
  //   V  #V by 3 list of rest pose vertex positions
  //   T  #T list of bone pose transformations
  //   Omega #V by #T*10 list of precomputated matrix values
  // Outputs:
  //   U  #V by 3 list of output vertex positions
  template <
    typename DerivedV,
    typename DerivedOmega,
    typename DerivedU>
  IGL_INLINE void direct_delta_mush(
    const Eigen::MatrixBase<DerivedV> & V,
    const std::vector<
      Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>
    > & T, /* should eventually be templated more generally than double */
    const Eigen::MatrixBase<DerivedOmega> & Omega,
    Eigen::PlainObjectBase<DerivedU> & U);

  // Precomputation
  //
  // Inputs:
  //   V  #V by 3 list of rest pose vertex positions
  //   F  #F by 3 list of triangle indices into rows of V
  //   W  #V by #Edges list of weights
  //   p  number of smoothing iterations
  //   lambda  rotation smoothing step size
  //   kappa   translation smoothness step size
  //   alpha   translation smoothness blending weight
  // Outputs:
  //   Omega  #V by #T*10 list of precomputated matrix values
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedW,
    typename DerivedOmega>
  IGL_INLINE void direct_delta_mush_precomputation(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedW> & W,
    const int p,
    const typename DerivedV::Scalar lambda,
    const typename DerivedV::Scalar kappa,
    const typename DerivedV::Scalar alpha,
    Eigen::PlainObjectBase<DerivedOmega> & Omega);
} // namespace igl

#ifndef IGL_STATIC_LIBRARY
#  include "direct_delta_mush.cpp"
#endif

#endif
