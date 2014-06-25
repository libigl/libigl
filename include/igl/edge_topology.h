// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EDGE_TOPOLOGY_H
#define IGL_EDGE_TOPOLOGY_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  // Initialize Edges and their topological relations
  
  // Output:
  // EV  : #Ex2, Stores the edge description as pair of indices to vertices
  // FE : #Fx3, Stores the Triangle-Edge relation
  // EF : #Ex2: Stores the Edge-Triangle relation
template <typename DerivedV, typename DerivedF>
  IGL_INLINE void edge_topology(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F, 
    Eigen::MatrixXi& EV, 
    Eigen::MatrixXi& FE, 
    Eigen::MatrixXi& EF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "edge_topology.cpp"
#endif

#endif
