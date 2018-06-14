// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Gavin Barill <gavinpcb@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/

#ifndef IGL_POINT_AREAS_AND_NORMALS_H
#define IGL_POINT_AREAS_AND_NORMALS_H
#include "../../igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
    // Given a 3D set of points P, each with a list of k-nearest-neighbours,
    // estimate the geodesic voronoi area associated with each point.
    //
    // The k nearest neighbours may be known from running igl::knn_octree on
    // the output data from igl::build_octree. We reccomend using a k value
    // between 15 and 20 inclusive for accurate area estimation.
    //
    // If given, O represents an input oritentaiton to each point. This could
    // be an initial guess at the normals, a vector to a camera position (for
    // scanned data), or ground truth normals if you're only interested in the
    // area estimation. It's used to ensure area estimation only occurs using
    // neighbors that are on the same side of the surface (ie for thin
    // sheets), as well as to solve the orientation ambiguity of normal
    // estimation.
    //
    // Inputs:
    //   P  #P by 3 list of point locations
    //   I  #P by k list of k-nearest-neighbor indices into P
    //   O  #P by 3 list of point orientation vectors (optional)
    // Outputs:
    //   A  #P list of estimated areas
    //   N  #P by 3 list of estimated normals
    template <typename DerivedP, typename DerivedI, typename DerivedO,
      typename DerivedA, typename DerivedN>
    IGL_INLINE void point_areas_and_normals(
                                        const Eigen::MatrixBase<DerivedP>& P,
                                        const Eigen::MatrixBase<DerivedI>& I,
                                        const Eigen::MatrixBase<DerivedO>& O,
                                        Eigen::PlainObjectBase<DerivedA> & A,
                                        Eigen::PlainObjectBase<DerivedN> & N);
      
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "point_areas_and_normals.cpp"
#endif

#endif

