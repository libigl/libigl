// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CGAL_REMESH_INTERSECTIONS_H
#define IGL_CGAL_REMESH_INTERSECTIONS_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include "CGAL_includes.hpp"

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
    // Remesh faces according to results of intersection detection and
    // construction (e.g. from `igl::cgal::intersect_other` or
    // `igl::cgal::SelfIntersectMesh`)
    //
    // Inputs:
    //   V  #V by 3 list of vertex positions
    //   F  #F by 3 list of triangle indices into V
    //   T  #F list of cgal triangles
    //   offending #offending map taking face indices into F to pairs of order
    //     of first finding and list of intersection objects from all
    //     intersections
    //   edge2faces  #edges <= #offending*3 to incident offending faces 
    // Outputs:
    //   VV  #VV by 3 list of vertex positions
    //   FF  #FF by 3 list of triangle indices into V
    //   IF  #intersecting face pairs by 2  list of intersecting face pairs,
    //     indexing F
    //   J  #FF list of indices into F denoting birth triangle
    //   IM  #VV list of indices into VV of unique vertices.
    //
    template <
      typename DerivedV,
      typename DerivedF,
      typename Kernel,
      typename DerivedVV,
      typename DerivedFF,
      typename DerivedJ,
      typename DerivedIM>
    IGL_INLINE void remesh_intersections(
      const Eigen::PlainObjectBase<DerivedV> & V,
      const Eigen::PlainObjectBase<DerivedF> & F,
      const std::vector<CGAL::Triangle_3<Kernel> > & T,
      const std::map<
        typename DerivedF::Index,
        std::pair<typename DerivedF::Index,
          std::vector<CGAL::Object> > > & offending,
      const std::map<
        std::pair<typename DerivedF::Index,typename DerivedF::Index>,
        std::vector<typename DerivedF::Index> > & edge2faces,
      Eigen::PlainObjectBase<DerivedVV> & VV,
      Eigen::PlainObjectBase<DerivedFF> & FF,
      Eigen::PlainObjectBase<DerivedJ> & J,
      Eigen::PlainObjectBase<DerivedIM> & IM);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "remesh_intersections.cpp"
#endif
  
#endif

