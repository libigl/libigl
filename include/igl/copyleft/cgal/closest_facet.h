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
      //
      // Outputs:
      //   R  #P list of closest facet indices.
      //   S  #P list of bools indicating on which side of the closest facet
      //      each query point lies.
      template<
          typename DerivedV,
          typename DerivedF,
          typename DerivedI,
          typename DerivedP,
          typename uE2EType,
          typename DerivedEMAP,
          typename DerivedR,
          typename DerivedS >
      IGL_INLINE void closest_facet(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedI>& I,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const std::vector<std::vector<uE2EType> >& uE2E,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
              Eigen::PlainObjectBase<DerivedR>& R,
              Eigen::PlainObjectBase<DerivedS>& S);

      template<
          typename DerivedV,
          typename DerivedF,
          typename DerivedP,
          typename uE2EType,
          typename DerivedEMAP,
          typename DerivedR,
          typename DerivedS >
      IGL_INLINE void closest_facet(
              const Eigen::PlainObjectBase<DerivedV>& V,
              const Eigen::PlainObjectBase<DerivedF>& F,
              const Eigen::PlainObjectBase<DerivedP>& P,
              const std::vector<std::vector<uE2EType> >& uE2E,
              const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
              Eigen::PlainObjectBase<DerivedR>& R,
              Eigen::PlainObjectBase<DerivedS>& S);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#include "closest_facet.cpp"
#endif
#endif
