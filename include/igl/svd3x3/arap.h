// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ARAP_H
#define IGL_ARAP_H
#include <igl/igl_inline.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/ARAPEnergyType.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
  struct ARAPData
  {
    // n  #V
    // G  #V list of group indices (1 to k) for each vertex, such that vertex i
    //    is assigned to group G(i)
    // energy  type of energy to use
    // with_dynamics  whether using dynamics 
    // f_ext  #V by dim list of external forces
    // h  dynamics time step
    // max_iter  maximum inner iterations
    // K  rhs pre-multiplier
    // solver_data  quadratic solver data
    // b  list of boundary indices into V
    int n;
    Eigen::VectorXi G;
    ARAPEnergyType energy;
    bool with_dynamics;
    Eigen::MatrixXd f_ext;
    double h;
    int max_iter;
    Eigen::SparseMatrix<double> K;
    Eigen::SparseMatrix<double> CSM;
    min_quad_with_fixed_data<double> solver_data;
    Eigen::VectorXi b;
      ARAPData():
        n(0),
        G(),
        energy(ARAP_ENERGY_TYPE_DEFAULT),
        with_dynamics(false),
        f_ext(),
        h(1),
        max_iter(10),
        K(),
        CSM(),
        solver_data(),
        b()
    {
    };
  };
  
  // Compute necessary information to start using an ARAP deformation
  //
  // Inputs:
  //   V  #V by dim list of mesh positions
  //   F  #F by simplex-size list of triangle|tet indices into V
  //   b  #b list of "boundary" fixed vertex indices into V
  // Outputs:
  //   data  struct containing necessary precomputation
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedb>
  IGL_INLINE bool arap_precomputation(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<Derivedb> & b,
    ARAPData & data);

  template <
    typename Derivedbc,
    typename DerivedU>
  IGL_INLINE bool arap_solve(
    const Eigen::PlainObjectBase<Derivedbc> & bc,
    ARAPData & data,
    Eigen::PlainObjectBase<DerivedU> & U);
};

#ifdef IGL_HEADER_ONLY
#include "arap.cpp"
#endif

#endif
