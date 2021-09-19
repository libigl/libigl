// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_REMESH_SELF_INTERSECTIONS_PARAM_H
#define IGL_COPYLEFT_CGAL_REMESH_SELF_INTERSECTIONS_PARAM_H

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Optional Parameters
      //   DetectOnly  Only compute IF, leave VV and FF alone
      //
      // detect_only  avoid constructing intersections results when possible
      // first_only  return after detecting the first intersection (if
      //   first_only==true, then detect_only should also be true)
      // stitch_all  whether to stitch all resulting constructed elements into a
      //   (non-manifold) mesh 
      // slow_and_more_precise_rounding  whether to use slow and more precise
      //   rounding (see assign_scalar)
      struct RemeshSelfIntersectionsParam
      {
        bool detect_only;
        bool first_only;
        bool stitch_all;
        bool slow_and_more_precise_rounding;
        inline RemeshSelfIntersectionsParam(
          bool _detect_only=false, 
          bool _first_only=false,
          bool _stitch_all=false,
          bool _slow_and_more_precise_rounding=false
          ):
          detect_only(_detect_only),
          first_only(_first_only),
          stitch_all(_stitch_all),
          slow_and_more_precise_rounding(_slow_and_more_precise_rounding)
        {};
      };
    }
  }
}

#endif
