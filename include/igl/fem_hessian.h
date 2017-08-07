// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com> and Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_FEM_HESSIAN_H
#define IGL_FEM_HESSIAN_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace igl
{
    // Constructs the finite element Hessian matrix
    // as described in https://arxiv.org/abs/1707.04348
    // The interior vertices are NOT set to zero yet.
    //
    // Inputs:
    //   V  #V by dim list of mesh vertex positions
    //   F  #F by 3 list of mesh faces (must be triangles)
    // Outputs:
    //   H  #V by #V Hessian energy matrix, each column i corresponding to V(i,:)
    //
    //
    //
    template <typename DerivedV, typename DerivedF, typename Scalar>
    IGL_INLINE void fem_hessian(
                                const Eigen::PlainObjectBase<DerivedV> & V,
                                const Eigen::PlainObjectBase<DerivedF> & F,
                                Eigen::SparseMatrix<Scalar>& H);
    
}

#ifndef IGL_STATIC_LIBRARY
#  include "fem_hessian.cpp"
#endif

#endif
