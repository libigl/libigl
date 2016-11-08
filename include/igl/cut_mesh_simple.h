// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Christian Schüller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CUT_MESH_SIMPLE
#define IGL_CUT_MESH_SIMPLE
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

#include <igl/HalfEdgeIterator.h>

namespace igl
{

  // Description
  // 
  //
  // Inputs:
  //   V:
  //   T:
  //   cut: needs to be a valid vertex sequence 
  //
  // Outputs:
  //
  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut);

  template <typename DerivedS, typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices);

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<int>& cutVerticesLink);

  template <typename DerivedS,typename DerivedI,int Option>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<int>& cutVerticesLink,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutHalfedges,
                std::vector<int>& cutHalfedgesLink);
  
  template <typename DerivedS,typename DerivedI,int Option>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<std::vector<int>>& cutVerticesLink);

  template <typename DerivedS,typename DerivedI,int Option>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<std::vector<int>>& cutVerticesLink,
                std::vector<Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,4,Option>>& cutHalfedges,
                std::vector<std::vector<int>>& cutHalfedgesLink);
};


#ifndef IGL_STATIC_LIBRARY
#include "cut_mesh_simple.cpp"
#endif


#endif
