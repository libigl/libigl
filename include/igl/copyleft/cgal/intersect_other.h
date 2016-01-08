// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_INTERSECT_OTHER_H
#define IGL_COPYLEFT_CGAL_INTERSECT_OTHER_H
#include "../../igl_inline.h"
#include "RemeshSelfIntersectionsParam.h"

#include <Eigen/Dense>

#ifdef MEX
#  include <mex.h>
#  include <cassert>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // INTERSECT_OTHER Given a triangle mesh (VA,FA) and another mesh (VB,FB)
      // find all pairs of intersecting faces. Note that self-intersections are
      // ignored.
      // 
      // Inputs:
      //   VA  #V by 3 list of vertex positions
      //   FA  #F by 3 list of triangle indices into VA
      //   VB  #V by 3 list of vertex positions
      //   FB  #F by 3 list of triangle indices into VB
      //   params   whether to detect only and then whether to only find first
      //     intersection
      // Outputs:
      //   IF  #intersecting face pairs by 2 list of intersecting face pairs,
      //     indexing FA and FB
      //   VVA  #VVA by 3 list of vertex positions
      //   FFA  #FFA by 3 list of triangle indices into VVA
      //   JA  #FFA list of indices into FA denoting birth triangle
      //   IMA  #VVA list of indices into VVA of unique vertices.
      //   VVB  #VVB by 3 list of vertex positions
      //   FFB  #FFB by 3 list of triangle indices into VVB
      //   JB  #FFB list of indices into FB denoting birth triangle
      //   IMB  #VVB list of indices into VVB of unique vertices.
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename DerivedIF,
        typename DerivedVVA,
        typename DerivedFFA,
        typename DerivedJA,
        typename DerivedIMA,
        typename DerivedVVB,
        typename DerivedFFB,
        typename DerivedJB,
        typename DerivedIMB>
      IGL_INLINE bool intersect_other(
        const Eigen::PlainObjectBase<DerivedVA> & VA,
        const Eigen::PlainObjectBase<DerivedFA> & FA,
        const Eigen::PlainObjectBase<DerivedVB> & VB,
        const Eigen::PlainObjectBase<DerivedFB> & FB,
        const RemeshSelfIntersectionsParam & params,
        Eigen::PlainObjectBase<DerivedIF> & IF,
        Eigen::PlainObjectBase<DerivedVVA> & VVA,
        Eigen::PlainObjectBase<DerivedFFA> & FFA,
        Eigen::PlainObjectBase<DerivedJA>  & JA,
        Eigen::PlainObjectBase<DerivedIMA> & IMA,
        Eigen::PlainObjectBase<DerivedVVB> & VVB,
        Eigen::PlainObjectBase<DerivedFFB> & FFB,
        Eigen::PlainObjectBase<DerivedJB>  & JB,
        Eigen::PlainObjectBase<DerivedIMB> & IMB);
      // Legacy wrapper for detect only using common types.
      //
      // Inputs:
      //   VA  #V by 3 list of vertex positions
      //   FA  #F by 3 list of triangle indices into VA
      //   VB  #V by 3 list of vertex positions
      //   FB  #F by 3 list of triangle indices into VB
      //   first_only  whether to only detect the first intersection.
      // Outputs:
      //   IF  #intersecting face pairs by 2 list of intersecting face pairs,
      //     indexing FA and FB
      // Returns true if any intersections were found
      //
      // See also: remesh_self_intersections
      IGL_INLINE bool intersect_other(
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        const Eigen::MatrixXi & FB,
        const bool first_only,
        Eigen::MatrixXi & IF);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "intersect_other.cpp"
#endif
  
#endif

