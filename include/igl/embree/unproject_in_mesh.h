// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_UNPROJECT_IN_MESH
#define IGL_UNPROJECT_IN_MESH
#include <igl/igl_inline.h>
#include <Eigen/Core>

#include <vector>
#include "Hit.h"

namespace igl
{
  // Forward define
  class EmbreeIntersector;

  #ifndef IGL_OPENGL_4
  // Unproject a screen location (using current opengl viewport, projection, and
  // model view) to a 3D position
  //
  // Inputs:
  //    x  x-coordinate of mouse location
  //    y  y-coordinate of mouse location
  //    ei  EmbreeIntersector containing (V,F)
  // Outputs:
  //    obj  3d unprojected mouse point in mesh
  // Returns number of hits
  //
  template <
    typename Derivedobj>
  IGL_INLINE int unproject_in_mesh(
    const double x,
    const double y,
    const igl::EmbreeIntersector & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj);

  template <
    typename Derivedobj>
  IGL_INLINE int unproject_in_mesh(
    const double x,
    const double y,
    const igl::EmbreeIntersector & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj,
    std::vector<igl::Hit > & hits);
  #endif

// Unproject a screen location (using the given model, proj and viewewport) to a 3D position
// and a set of hits
//
// Inputs:
//    pos        screen space coordinates
//    model      model matrix
//    proj       projection matrix
//    viewport   vieweport vector
//    ei         EmbreeIntersector containing (V,F)
// Outputs:
//    obj        3d unprojected mouse point in mesh
//    hits       vector of embree hits
// Returns number of hits
template <
  typename Derivedobj>
IGL_INLINE int unproject_in_mesh(
  const Eigen::Vector2f& pos,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj,
  std::vector<igl::Hit > & hits);

}
#ifndef IGL_STATIC_LIBRARY
#  include "unproject_in_mesh.cpp"
#endif
#endif
