// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_HAUSDORFF_H
#define IGL_HAUSDORFF_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl 
{
  // HAUSDORFF compute lower and upper bounds on the **directed** (asymmetric)
  // Hausdorff distance from mesh (VA,FA) to mesh (VB,FB) within a given error
  // bound epsilon:
  //
  // l(A,B) ≤ max min d(a,b) ≤ u(A,B)
  //          a∈A b∈B
  // u(A,B)-l(A,B) ≤ epsilon
  //
  // This is mostly an implementation of "Fast and accurate Hausdorff distance
  // calculation between meshes" [Guthe et al. 2005], with some ideas from
  // "Interactive Hausdorff Distance Computation for General Polygonal Models"
  // [Tang et al. 2009].
  // 
  // Inputs:
  //   VA  #VA by 3 list of vertex positions
  //   FA  #FA by 3 list of face indices into VA
  //   VB  #VB by 3 list of vertex positions
  //   FB  #FB by 3 list of face indices into VB
  //   eps  desired bound on error
  // Outputs:
  //   lower  lower bound on directed Hausdorff distance
  //   upper  upper bound on directed Hausdorff distance
  template <
    typename DerivedVA, 
    typename DerivedFA,
    typename DerivedVB,
    typename DerivedFB,
    typename Scalar>
  IGL_INLINE void hausdorff(
    const Eigen::PlainObjectBase<DerivedVA> & OVA, 
    const Eigen::PlainObjectBase<DerivedFA> & OFA,
    const Eigen::PlainObjectBase<DerivedVB> & VB, 
    const Eigen::PlainObjectBase<DerivedFB> & FB,
    const Scalar eps,
    Scalar & lower,
    Scalar & upper);
  // HAUSDORFF compute the **symmetric** Hausdorff distance between mesh
  // (VA,FA) and mesh (VB,FB). That is:
  //
  // d(A,B) = max ( max min d(a,b) , max min d(b,a) )
  //                a∈A b∈B          b∈B a∈A
  //
  // Note: this is an iterative method that will compute the distance within
  // double precision epsilon times the bounding box diagonal of the combined
  // meshes A and B.
  //
  // Inputs:
  //   VA  #VA by 3 list of vertex positions
  //   FA  #FA by 3 list of face indices into VA
  //   VB  #VB by 3 list of vertex positions
  //   FB  #FB by 3 list of face indices into VB
  // Outputs:
  //   d  hausdorff distance
  //
  template <
    typename DerivedVA, 
    typename DerivedFA,
    typename DerivedVB,
    typename DerivedFB,
    typename Scalar>
  IGL_INLINE void hausdorff(
    const Eigen::PlainObjectBase<DerivedVA> & VA, 
    const Eigen::PlainObjectBase<DerivedFA> & FA,
    const Eigen::PlainObjectBase<DerivedVB> & VB, 
    const Eigen::PlainObjectBase<DerivedFB> & FB,
    Scalar & d);
}

#ifndef IGL_STATIC_LIBRARY
#  include "hausdorff.cpp"
#endif

#endif

