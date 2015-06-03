// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_POLYVECTOR_FIELD_COMB_FROM_MATCHINGS_AND_CUTS
#define IGL_POLYVECTOR_FIELD_COMB_FROM_MATCHINGS_AND_CUTS
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl {
  // Given an mesh, an N- polyvector field with given matchings between
  // the vector sets across interior edges, and corresponding mesh cuts,
  // compute a "combed" polyvector field, i.e.
  // separate the vector set field into N vector fields, where the separation is defined
  // by the matchings. The cuts have previously been generated based on the field
  // singularities, see igl::polyvector_field_cut_mesh_with_singularities.
  // Inputs:
  //   V                #V by 3 list of the vertex positions
  //   F                #F by 3 list of the faces (must be triangles)
  //   TT               #F by 3 triangle to triangle adjacent matrix (e.g. computed
  //                    via igl:triangle_triangle_adjacency).
  //   E2F              #E by 2 list of the edge-to-face relation (e.g. computed
  //                    via igl::edge_topology)
  //   F2E              #F by 3 list of the face-to-edge relation (e.g. computed
  //                    via igl::edge_topology)
  //   sol3D            #F by 3n list of the 3D coordinates of the per-face vectors of the input vector set field
  //                    (stacked horizontally for each triangle). Vector #1 in one face does not necessarily match
  //                    vector #1 in the adjacent face.
  //   match_ab         #E by N matrix, describing for each edge the matching a->b, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #i of
  //                    the vector set in a is matched to vector #mab[i] in b)
  //   match_ba         #E by N matrix, describing for each edge the matching b->a, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #mba[i] of
  //                    the vector set in a is matched to vector #i in b)
  //   cuts             #F by 3 list of boolean flags, indicating the edges that need to be cut
  // Outputs:
  //   sol3D_combed     #F by 3n list of the 3D coordinates of the per-face vectors of the combed vector set field
  //                    (stacked horizontally for each triangle). Vector #1 in one face will match vector #1 in
  //                    the adjacent face.
  //
  template <typename DerivedV, typename DerivedF, typename DerivedTT, typename DerivedS, typename DerivedM, typename DerivedC>
  IGL_INLINE void polyvector_field_comb_from_matchings_and_cuts(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const Eigen::PlainObjectBase<DerivedTT> &TT,
                                                                const Eigen::MatrixXi &E2F,
                                                                const Eigen::MatrixXi &F2E,
                                                                const Eigen::PlainObjectBase<DerivedS> &sol3D,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                const Eigen::PlainObjectBase<DerivedC> &cuts,
                                                                Eigen::PlainObjectBase<DerivedS> &sol3D_combed);
  
  //Wrapper with only the mesh as input
  template <typename DerivedV, typename DerivedF, typename DerivedS, typename DerivedM, typename DerivedC>
  IGL_INLINE void polyvector_field_comb_from_matchings_and_cuts(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const Eigen::PlainObjectBase<DerivedS> &sol3D,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                const Eigen::PlainObjectBase<DerivedC> &cuts,
                                                                Eigen::PlainObjectBase<DerivedS> &sol3D_combed);
};


#ifndef IGL_STATIC_LIBRARY
#include "polyvector_field_comb_from_matchings_and_cuts.cpp"
#endif


#endif /* defined(IGL_POLYVECTOR_FIELD_COMB_FROM_MATCHINGS_AND_CUTS) */
