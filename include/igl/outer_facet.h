// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OUTER_FACET_H
#define IGL_OUTER_FACET_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute a single facet which is guaranteed to be part of the "outer hull of
  // a mesh (V,F). This implementation follows Section 3.6 of "Direct repair of
  // self-intersecting meshes" [Attene 2014].
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  //   N  #F by 3 list of face normals
  //   I  #I list of facets to actually consider
  // Outputs:
  //   f  index of facet into V
  //   flip  whether facet's orientation should be flipped so that
  //     counter-clockwise normal points outward.
  //
  // See also: cgal/outer_hull.h
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename DerivedI,
    typename f_type>
  IGL_INLINE void outer_facet(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedN> & N,
    const Eigen::PlainObjectBase<DerivedI> & I,
    f_type & f,
    bool & flip);
}
#ifndef IGL_STATIC_LIBRARY
#  include "outer_facet.cpp"
#endif
#endif
