// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPEUP_LOCAL_PROJECTIONS_H
#define IGL_SHAPEUP_LOCAL_PROJECTIONS_H

#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <igl/cat.h>
#include <Eigen/Core>
#include <vector>


//This file implements several basic lcaol projection functions for the shapeup algorithm in shapeup.h

namespace igl
{
  
  //This projection does nothing but render points into projP. Mostly used for "echoing" the global step
  IGL_INLINE bool shapeup_identity_projection(const Eigen::PlainObjectBase<Eigen::MatrixXd>& P, const Eigen::PlainObjectBase<Eigen::VectorXi>& SC, const Eigen::PlainObjectBase<Eigen::MatrixXi>& S,  Eigen::PlainObjectBase<Eigen::MatrixXd>& projP);
  
  //the projection assumes that the sets are vertices of polygons in order
  IGL_INLINE bool shapeup_regular_face_projection(const Eigen::PlainObjectBase<Eigen::MatrixXd>& P, const Eigen::PlainObjectBase<Eigen::VectorXi>& SC, const Eigen::PlainObjectBase<Eigen::MatrixXi>& S,  Eigen::PlainObjectBase<Eigen::MatrixXd>& projP);
  
}

#ifndef IGL_STATIC_LIBRARY
#include "shapeup_local_projections.cpp"
#endif

#endif
