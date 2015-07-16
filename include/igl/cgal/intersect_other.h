// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CGAL_INTERSECT_OTHER_H
#define IGL_CGAL_INTERSECT_OTHER_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>

#ifdef MEX
#  include <mex.h>
#  include <cassert>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

namespace igl
{
  namespace cgal
  {
    // INTERSECT Given a triangle mesh (V,F) and another mesh (U,G) find all pairs
    // of intersecting faces. Note that self-intersections are ignored.
    // 
    // [VV,FF,IF] = selfintersect(V,F,'ParameterName',ParameterValue, ...)
    //
    // Inputs:
    //   V  #V by 3 list of vertex positions
    //   F  #F by 3 list of triangle indices into V
    //   U  #U by 3 list of vertex positions
    //   G  #G by 3 list of triangle indices into U
    //   first_only  whether to only detect the first intersection.
    // Outputs:
    //   IF  #intersecting face pairs by 2 list of intersecting face pairs,
    //     indexing F and G
    //
    // See also: selfintersect
    IGL_INLINE void intersect_other(
      const Eigen::MatrixXd & V,
      const Eigen::MatrixXi & F,
      const Eigen::MatrixXd & U,
      const Eigen::MatrixXi & G,
      const bool first_only,
      Eigen::MatrixXi & IF);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "intersect_other.cpp"
#endif
  
#endif

