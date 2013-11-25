// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "grad.h"

template <typename T, typename S>
IGL_INLINE void igl::grad(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
  const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> &F,
  const Eigen::Matrix<T, Eigen::Dynamic, 1>&X,
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &G )
{
  G = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(F.rows(),3);
  for (int i = 0; i<F.rows(); ++i)
  {
    // renaming indices of vertices of triangles for convenience
    int i1 = F(i,0);
    int i2 = F(i,1);
    int i3 = F(i,2);
    
    // #F x 3 matrices of triangle edge vectors, named after opposite vertices
    Eigen::Matrix<T, 1, 3> v32 = V.row(i3) - V.row(i2);
    Eigen::Matrix<T, 1, 3> v13 = V.row(i1) - V.row(i3);
    Eigen::Matrix<T, 1, 3> v21 = V.row(i2) - V.row(i1);
    
    // area of parallelogram is twice area of triangle
    // area of parallelogram is || v1 x v2 || 
    Eigen::Matrix<T, 1, 3> n  = v32.cross(v13); 
    
    // This does correct l2 norm of rows, so that it contains #F list of twice
    // triangle areas
    double dblA = std::sqrt(n.dot(n));
    
    // now normalize normals to get unit normals
    Eigen::Matrix<T, 1, 3> u = n / dblA;
    
    // rotate each vector 90 degrees around normal
    double norm21 = std::sqrt(v21.dot(v21));
    double norm13 = std::sqrt(v13.dot(v13));
    Eigen::Matrix<T, 1, 3> eperp21 = u.cross(v21);
    eperp21 = eperp21 / std::sqrt(eperp21.dot(eperp21));
    eperp21 *= norm21;
    Eigen::Matrix<T, 1, 3> eperp13 = u.cross(v13);
    eperp13 = eperp13 / std::sqrt(eperp13.dot(eperp13));
    eperp13 *= norm13;
    
    G.row(i) = ((X[i2] -X[i1]) *eperp13 + (X[i3] - X[i1]) *eperp21) / dblA;
  };
}
  
  

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
