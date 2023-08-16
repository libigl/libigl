// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COMISO_NROSY_H
#define IGL_COMISO_NROSY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../../igl_inline.h"

namespace igl
{
  namespace copyleft
  {
  namespace comiso
  {
    /// Generate a N-RoSy field from a sparse set of constraints
    ///
    /// @param[in] V       #V by 3 list of mesh vertex coordinates
    /// @param[in] F       #F by 3 list of mesh faces (must be triangles)
    /// @param[in] b       #B by 1 list of constrained face indices
    /// @param[in] bc      #B by 3 list of representative vectors for the
    ///   constrained faces
    /// @param[in] b_soft  #S by 1 b for soft constraints
    /// @param[in] w_soft  #S by 1 weight for the soft constraints (0-1)
    /// @param[in] bc_soft #S by 3 bc for soft constraints
    /// @param[in] N       the degree of the N-RoSy vector field
    /// @param[in] soft    the strength of the soft constraints w.r.t. smoothness
    ///           (0 -> smoothness only, 1->constraints only)
    /// @param[out] R       #F by 3 the representative vectors of the interpolated field
    /// @param[out] S       #V by 1 the singularity index for each vertex (0 = regular)
    IGL_INLINE void nrosy(
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      const Eigen::VectorXi& b,
      const Eigen::MatrixXd& bc,
      const Eigen::VectorXi& b_soft,
      const Eigen::VectorXd& w_soft,
      const Eigen::MatrixXd& bc_soft,
      int N,
      double soft,
      Eigen::MatrixXd& R,
      Eigen::VectorXd& S
      );
    /// \overload
    IGL_INLINE void nrosy(
     const Eigen::MatrixXd& V,
     const Eigen::MatrixXi& F,
     const Eigen::VectorXi& b,
     const Eigen::MatrixXd& bc,
     int N,
     Eigen::MatrixXd& R,
     Eigen::VectorXd& S
      );

  }
}
}

#ifndef IGL_STATIC_LIBRARY
#  include "nrosy.cpp"
#endif

#endif
