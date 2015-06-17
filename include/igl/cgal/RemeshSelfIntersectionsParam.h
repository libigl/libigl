// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CGAL_REMESH_SELF_INTERSECTIONS_PARAM_H
#define IGL_CGAL_REMESH_SELF_INTERSECTIONS_PARAM_H

namespace igl
{
  namespace cgal
  {
    // Optional Parameters
    //   DetectOnly  Only compute IF, leave VV and FF alone
    struct RemeshSelfIntersectionsParam
    {
      bool detect_only;
      bool first_only;
      RemeshSelfIntersectionsParam():detect_only(false),first_only(false){};
    };
  }
}

#endif
