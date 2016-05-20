// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Amir Vaxman <avaxman@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_GENERAL_EDGE_TOPOLOGY_H
#define IGL_GENERAL_EDGE_TOPOLOGY_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  // Initialize Edges and their topological relations
    
    //input:
    //gD: #Fx1, Stores the degrees of the faces
    //gF: #Fx(max(D)), Stores the faces
  
  // Output:
  // EV  : #Ex2, Stores the edge description as pair of indices to vertices
  // FE : #Fx3, Stores the Face-Edge relation
  // EF : #Ex2: Stores the Edge-Face relation
template <typename DerivedF>
  IGL_INLINE void general_edge_topology(
    const Eigen::VectorXi& gD,
    const Eigen::PlainObjectBase<DerivedF>& gF,
    Eigen::MatrixXi& gEV,
    Eigen::MatrixXi& gFE,
    Eigen::MatrixXi& gEF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "general_edge_topology.cpp"
#endif

#endif
