// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PROCRUSTES_H
#define IGL_PROCRUSTES_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace igl
{
    // Solve Procrustes problem in d dimensions.
    // Given two point sets X,Y in R^d find best scale s, rotation/reflection R  and translation t 
    // s.t. |s*X*R + t - Y|^2 is minimized.
    //
    // Example:
    // MatrixXd X, Y; (containing 3d points as rows)
    // AffineCompact3d T;
    // igl::procrustes(X,Y,true,false,T);
    // MatrixXd Xprime = (X * T.linear()).rowwise() + T.translation().transpose();
    //     
    //
    // Templates:
    //    DerivedV point type
    //    Scalar   scalar type
    //    DIM      point dimension
    //    TType    type of transformation (Isometry,Affine,AffineCompact,Projective)
    // Inputs:
    //    X  #V by DIM first list of points
    //    Y  #V by DIM second list of points
    //    includeScaling  if scaling should be allowed
    //    includeReflections  if R is allowed to be a reflection
    // Outputs:
    //    T  transformation that minimizes error    
    template <typename DerivedV, typename Scalar, int DIM, int TType>
    IGL_INLINE void procrustes(
        const Eigen::PlainObjectBase<DerivedV>& X,
        const Eigen::PlainObjectBase<DerivedV>& Y,
        bool includeScaling,
        bool includeReflections,
        Eigen::Transform<Scalar,DIM,TType>& T);
}

#ifndef IGL_STATIC_LIBRARY
    #include "procrustes.cpp"
#endif

#endif