// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_TRIANGLE_CLOSEST_FEATURE_H
#define IGL_PREDICATES_EXACT_POINT_TRIANGLE_CLOSEST_FEATURE_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  enum class ClosestPointFeature { VERTEX=0, EDGE=1, FACE=2 };
  /// Classify the unique closest point on a closed triangle.
  /// Feature ids are corner i, open edge [i,(i+1)%3], or 0 for the open face.
  /// This edge convention matches the array triangle_triangle_adjacency API.
  /// Degenerate triangles are supported: prefer the first coincident vertex,
  /// otherwise the first open edge containing the closest point.
  /// Inputs must be finite float or double coordinates. Link igl::predicates
  /// with its double backend. Correctness assumes round-to-nearest with gradual
  /// underflow, without reassociation or contraction. Exact refinement rescales
  /// coordinates by a common power of two and throws std::range_error if their
  /// exponent span cannot fit safely. A sufficient range for the entire point
  /// predicate family is: after rescaling, every nonzero coordinate's lowest
  /// set bit is at least 2^-100 and its magnitude is less than 2^91. Predicates
  /// of lower degree can refine wider spans. No expansion survives the call.
  /// @param[in] query,a,b,c  Query and triangle corner 3-vectors
  /// @param[out] feature_id  Local corner, edge, or face id as above
  /// @return Dimension of the feature containing the closest point
  template<typename DerivedQ,typename DerivedA>
  IGL_INLINE ClosestPointFeature exact_point_triangle_closest_feature(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedA> &a,
    const Eigen::MatrixBase<DerivedA> &b,
    const Eigen::MatrixBase<DerivedA> &c,int &feature_id);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "exact_point_triangle_closest_feature.cpp"
#endif
#endif
