// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BIJECTIVE_COMPOSITE_HARMONIC_MAPPING_H
#define IGL_BIJECTIVE_COMPOSITE_HARMONIC_MAPPING_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl 
{
  // Compute a planar mapping of a triangulated polygon (V,F) subjected to
  // boundary conditions (b,bc). The mapping should be bijective in the sense
  // that no triangles' areas become negative (this assumes they started
  // positive). This mapping is computed by "composing" harmonic mappings
  // between incremental morphs of the boundary conditions. This is a bit like
  // a discrete version of "Bijective Composite Mean Value Mappings" [Schneider
  // et al. 2013] but with a discrete harmonic map (cf. harmonic coordinates)
  // instead of mean value coordinates. This is inspired by "Embedding a
  // triangular graph within a given boundary" [Xu et al. 2011].
  //
  // Inputs:
  //   V  #V by 2 list of triangle mesh vertex positions
  //   F  #F by 3 list of triangle indices into V
  //   b  #b list of boundary indices into V
  //   bc  #b by 2 list of boundary conditions corresponding to b
  // Outputs:
  //   U  #V by 2 list of output mesh vertex locations
  // Returns true if and only if U contains a successfull bijectie mapping
  //
  // 
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedb,
    typename Derivedbc,
    typename DerivedU>
  IGL_INLINE bool bijective_composite_harmonic_mapping(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<Derivedb> & b,
    const Eigen::MatrixBase<Derivedbc> & bc,
    Eigen::PlainObjectBase<DerivedU> & U);
}

#ifndef IGL_STATIC_LIBRARY
#  include "bijective_composite_harmonic_mapping.cpp"
#endif
#endif
