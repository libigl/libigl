// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_UNPROJECT_ONTO_MESH
#define IGL_UNPROJECT_ONTO_MESH
#include <igl/igl_inline.h>
#include <Eigen/Core>

#include <vector>
#include "Hit.h"

namespace igl
{
  // Forward define
  class EmbreeIntersector;
  // Unproject a screen location (using the given model, proj and viewport) to find
  // the first face on the mesh and the closest vertex
  //
  // Inputs:
  //    pos        screen space coordinates
  //    F          #F by 3 face matrix
  //    model      model matrix
  //    proj       projection matrix
  //    viewport   vieweport vector
  //    ei         EmbreeIntersector containing (V,F)
  // Outputs:
  //    fid        id of the first face hit
  //    vid        vertex id of the closest vertex hit
  // Returns true if there is a hit
  //
  // Known bugs: This should return the barycentric coordinates in F(fid,:)
  // rather than vid.
  IGL_INLINE bool unproject_onto_mesh(
    const Eigen::Vector2f& pos,
    const Eigen::MatrixXi& F,
    const Eigen::Matrix4f& model,
    const Eigen::Matrix4f& proj,
    const Eigen::Vector4f& viewport,
    const igl::EmbreeIntersector & ei,
    int& fid,
    int& vid);
}
#ifndef IGL_STATIC_LIBRARY
#  include "unproject_onto_mesh.cpp"
#endif
#endif
