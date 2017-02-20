// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Christian Schüller <schuellchr@gmail.com>
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

  // Description:
  // Cuts a triangle mesh along edges. 
 
  // Inputs:
  //   V: #V by 3 list of the vertex positions
  //   T: #F by 3 list of the faces (must be triangles)
  //   cut: A sequence of vertex ids defining the cut along edges (e.g. {1,2,3} defines two edges: {1,2} and {2,3})
  //   cuts: A list with multiple cuts. Provided edges don't need to be unique, can intersect each other or even be boundary edges.
  //   
  // Outputs:
  //  cutVertices: list of correspponding unique vertex id pairs after cutting
  //  cutVerticesLink: indices of 'cutVertices' corresponding to vertices in 'cut' 
  //  cutHalfedges: list of correspponding unique halfedges [triangleAId,edgeInAId,triangleBId,edgeInBId] (e.g. T(triangleAId,edgeInAId))
  //  cutHalfedgesLink: indices of 'cutHalfedges' corresponding to vertices in 'cut' (1. vertex id in 'cut' corresponds to 1. halfedge in 'cutHalfedges')
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
