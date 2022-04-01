// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef IGL_COPYLET_CGAL_CLOSEST_FACET_H
#define IGL_COPYLET_CGAL_CLOSEST_FACET_H

#include "../../igl_inline.h"
#include <Eigen/Core>
#include <vector>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Determine the closest facet for each of the input points.
      //
      // Inputs:
      //   V  #V by 3 array of vertices.
      //   F  #F by 3 array of faces.
      //   I  #I list of triangle indices to consider.
      //   P  #P by 3 array of query points.
      //  EMAP  #F*3 list of indices into uE.
      //  uEC  #uE+1 list of cumsums of directed edges sharing each unique edge
      //  uEE  #E list of indices into E (see `igl::unique_edge_map`)
      //   VF  #V list of lists of incident faces (adjacency list)
      //   VFi  #V list of lists of index of incidence within incident faces
      //     listed in VF
      //   tree  AABB containing triangles of (V,F(I,:))
      //   triangles  #I list of cgal triangles
      //   in_I  #F list of whether in submesh
      // Outputs:
      //   R  #P list of closest facet indices.
      //   S  #P list of bools indicating on which side of the closest facet
      //      each query point lies.
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename DerivedP,
        typename DerivedEMAP,
        typename DeriveduEC,
        typename DeriveduEE,
        typename Kernel,
        typename DerivedR,
        typename DerivedS >
      IGL_INLINE void closest_facet(
          const Eigen::PlainObjectBase<DerivedV>& V,
          const Eigen::PlainObjectBase<DerivedF>& F,
          const Eigen::PlainObjectBase<DerivedI>& I,
          const Eigen::PlainObjectBase<DerivedP>& P,
          const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
          const Eigen::PlainObjectBase<DeriveduEC>& uEC,
          const Eigen::PlainObjectBase<DeriveduEE>& uEE,
          const std::vector<std::vector<size_t> > & VF,
          const std::vector<std::vector<size_t> > & VFi,
          const CGAL::AABB_tree<
            CGAL::AABB_traits<
              Kernel, 
              CGAL::AABB_triangle_primitive<
                Kernel, typename std::vector<
                  typename Kernel::Triangle_3 >::iterator > > > & tree,
          const std::vector<typename Kernel::Triangle_3 > & triangles,
          const std::vector<bool> & in_I,
          Eigen::PlainObjectBase<DerivedR>& R,
          Eigen::PlainObjectBase<DerivedS>& S);
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename DerivedP,
        typename DerivedEMAP,
        typename DeriveduEC,
        typename DeriveduEE,
        typename DerivedR,
        typename DerivedS >
      IGL_INLINE void closest_facet(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedI>& I,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const Eigen::PlainObjectBase<DeriveduEC>& uEC,
        const Eigen::PlainObjectBase<DeriveduEE>& uEE,
              Eigen::PlainObjectBase<DerivedR>& R,
              Eigen::PlainObjectBase<DerivedS>& S);
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedP,
        typename DerivedEMAP,
        typename DeriveduEC,
        typename DeriveduEE,
        typename DerivedR,
        typename DerivedS >
      IGL_INLINE void closest_facet(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const Eigen::PlainObjectBase<DeriveduEC>& uEC,
        const Eigen::PlainObjectBase<DeriveduEE>& uEE,
        Eigen::PlainObjectBase<DerivedR>& R,
        Eigen::PlainObjectBase<DerivedS>& S);
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename DerivedP,
        typename DerivedEMAP,
        typename DeriveduEC,
        typename DeriveduEE,
        typename Kernel,
        typename DerivedR,
        typename DerivedS >
      IGL_INLINE void closest_facet(
          const Eigen::PlainObjectBase<DerivedV>& V,
          const Eigen::PlainObjectBase<DerivedF>& F,
          const Eigen::PlainObjectBase<DerivedI>& I,
          const Eigen::PlainObjectBase<DerivedP>& P,
          const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
          const Eigen::PlainObjectBase<DeriveduEC>& uEC,
          const Eigen::PlainObjectBase<DeriveduEE>& uEE,
          const std::vector<std::vector<size_t> > & VF,
          const std::vector<std::vector<size_t> > & VFi,
          const CGAL::AABB_tree<
            CGAL::AABB_traits<
              Kernel, 
              CGAL::AABB_triangle_primitive<
                Kernel, typename std::vector<
                  typename Kernel::Triangle_3 >::iterator > > > & tree,
          const std::vector<typename Kernel::Triangle_3 > & triangles,
          const std::vector<bool> & in_I,
          Eigen::PlainObjectBase<DerivedR>& R,
          Eigen::PlainObjectBase<DerivedS>& S);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#include "closest_facet.cpp"
#endif
#endif
