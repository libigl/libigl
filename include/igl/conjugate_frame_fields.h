// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_CONJUGATE_FRAME_FIELDS_H
#define IGL_CONJUGATE_FRAME_FIELDS_H
#include "igl_inline.h"
#include "ConjugateFFSolverData.h"
#include <Eigen/Core>
#include <vector>

namespace igl {
  //todo
  // TODO: isConstrained should become a list of indices for consistency with
  //       n_polyvector
  
  template <typename DerivedV, typename DerivedF, typename DerivedO>
  IGL_INLINE void conjugate_frame_fields(const Eigen::PlainObjectBase<DerivedV> &V,
                                         const Eigen::PlainObjectBase<DerivedF> &F,
                                         const Eigen::VectorXi &isConstrained,
                                         const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                         Eigen::PlainObjectBase<DerivedO> &output,
                                         int _maxIter = 50,
                                         const typename DerivedV::Scalar _lambdaOrtho = .1,
                                         const typename DerivedV::Scalar _lambdaInit = 100,
                                         const typename DerivedV::Scalar _lambdaMultFactor = 1.01,
                                         bool _doHardConstraints = true);
  
  template <typename DerivedV, typename DerivedF, typename DerivedO>
  IGL_INLINE typename DerivedV::Scalar conjugate_frame_fields(const ConjugateFFSolverData<DerivedV, DerivedF> &csdata,
                                                              const Eigen::VectorXi &isConstrained,
                                                              const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                                              Eigen::PlainObjectBase<DerivedO> &output,
                                                              int _maxIter = 50,
                                                              const typename DerivedV::Scalar _lambdaOrtho = .1,
                                                              const typename DerivedV::Scalar _lambdaInit = 100,
                                                              const typename DerivedV::Scalar _lambdaMultFactor = 1.01,
                                                              bool _doHardConstraints = true);
  
};


#ifndef IGL_STATIC_LIBRARY
#include "conjugate_frame_fields.cpp"
#endif


#endif
