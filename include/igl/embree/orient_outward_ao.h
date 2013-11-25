// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ORIENT_OUTWARD_AO_H
#define IGL_ORIENT_OUTWARD_AO_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Orient each component (identified by C) of a mesh (V,F) using ambient occlusion 
  // such that the front side is less occluded than back side
  //
  // Inputs:
  //   V                            #V by 3 list of vertex positions
  //   F                            #F by 3 list of triangle indices
  //   C                            #F list of components
  //   min_num_rays_per_component   Each component receives at least this number of rays
  //   total_num_rays               Total number of rays that will be shot
  // Outputs:
  //   FF  #F by 3 list of new triangle indices such that FF(~I,:) = F(~I,:) and
  //     FF(I,:) = fliplr(F(I,:)) (OK if &FF = &F)
  //   I  max(C)+1 list of whether face has been flipped
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedC, 
    typename DerivedFF, 
    typename DerivedI>
  IGL_INLINE void orient_outward_ao(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedC> & C,
    const int min_num_rays_per_component,
    const int total_num_rays,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
  
  // Call with default number of rays
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedC, 
    typename DerivedFF, 
    typename DerivedI>
  IGL_INLINE void orient_outward_ao(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
};

#ifdef IGL_HEADER_ONLY
#  include "orient_outward_ao.cpp"
#endif

#endif
