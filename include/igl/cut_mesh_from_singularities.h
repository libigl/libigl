// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_CUT_MESH_FROM_SINGULARITIES_H
#define IGL_CUT_MESH_FROM_SINGULARITIES_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  //todo
  
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedS, typename DerivedO>
  IGL_INLINE void cut_mesh_from_singularities(const Eigen::PlainObjectBase<DerivedV> &V,
                                                   const Eigen::PlainObjectBase<DerivedF> &F,
                                                   const Eigen::PlainObjectBase<DerivedM> &Handle_MMatch,
                                                   const Eigen::PlainObjectBase<DerivedS> &isSingularity,
                                                   const Eigen::PlainObjectBase<DerivedS> &singularityIndex,
                                                   Eigen::PlainObjectBase<DerivedO> &Handle_Seams);
}
#ifdef IGL_HEADER_ONLY
#include "cut_mesh_from_singularities.cpp"
#endif

#endif
