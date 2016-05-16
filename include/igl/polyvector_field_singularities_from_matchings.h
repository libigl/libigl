// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_POLYVECTOR_FIELD_SINGULARITIES_FROM_MATCHINGS
#define IGL_POLYVECTOR_FIELD_SINGULARITIES_FROM_MATCHINGS
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl {
  
  // We are given a polyvector field on a mesh (with N vectors per face), and matchings between the vector sets
  // across all non-boundary edges.  The function computes, for the one-ring
  // neighborhood of a given vertex, and for each of the vectors of the vector set,
  // the sequence of the vectors that match to it around the one-ring. If the vector that
  // we land on by following the matchings is not the original vector that we started from,
  // the vertex is a singularity.
  //
  // Inputs:
  //   V                #V by 3 list of the vertex positions
  //   F                #F by 3 list of the faces (must be triangles)
  //   VF               #V list of lists of incident faces (adjacency list), e.g.
  //                    as returned by igl::vertex_triangle_adjacency
  //   E2F              #E by 2 list of the edge-to-face relation (e.g. computed
  //                    via igl::edge_topology)
  //   F2E              #F by 3 list of the face-to-edge relation (e.g. computed
  //                    via igl::edge_topology)
  //   TT               #F by 3 triangle to triangle adjacent matrix (e.g. computed
  //                    via igl:triangle_triangle_adjacency)
  //   match_ab         #E by N matrix, describing for each edge the matching a->b, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #i of
  //                    the vector set in a is matched to vector #mab[i] in b)
  //   match_ba         #E by N matrix, describing for each edge the matching b->a, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #mba[i] of
  //                    the vector set in a is matched to vector #i in b)
  //   vi               the selected one ring
  //
  // Output:
  //   mvi              #numOneRingFaces by 1 list of the indices of the sequentially matching vectors
  //                    in the faces of the one ring (first enty is always vector_id, then the vector matching
  //                    vector_id in the next face, then the vector matching that in the third face etc.)
  //   fi               #numOneRingFaces by 1 list of the sequentially visited faces in the one ring neighborhood.
  //                    The one-ring is traversed in CLOCKWISE order with respect to the outward normal. (=opposite)
  //
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedTT>
  void polyvector_field_one_ring_matchings(const Eigen::PlainObjectBase<DerivedV> &V,
                                           const Eigen::PlainObjectBase<DerivedF> &F,
                                           const std::vector<std::vector<VFType> >& VF,
                                           const Eigen::MatrixXi& E2F,
                                           const Eigen::MatrixXi& F2E,
                                           const Eigen::PlainObjectBase<DerivedTT>& TT,
                                           const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                           const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                           const int vi,
                                           Eigen::MatrixXi &mvi,
                                           Eigen::VectorXi &fi);
  
  
  // Given a polyvector field on a mesh, and matchings between the vector sets
  // across all non-boundary edges, the function computes the singularities of the
  // polyvector field, by computing the one ring matchings.
  //
  // Inputs:
  //   V                #V by 3 list of the vertex positions
  //   F                #F by 3 list of the faces (must be triangles)
  //   V_border         #V by 1 list of booleans, indicating if the corresponging vertex is
  //                    at the mesh boundary, e.g. as returned by igl::is_border_vertex
  //   VF               #V list of lists of incident faces (adjacency list), e.g.
  //                    as returned by igl::vertex_triangle_adjacency
  //   TT               #F by 3 triangle to triangle adjacent matrix (e.g. computed
  //                    via igl:triangle_triangle_adjacency)
  //   E2F              #E by 2 list of the edge-to-face relation (e.g. computed
  //                    via igl::edge_topology)
  //   F2E              #F by 3 list of the face-to-edge relation (e.g. computed
  //                    via igl::edge_topology)
  //   match_ab         #E by N matrix, describing for each edge the matching a->b, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #i of
  //                    the vector set in a is matched to vector #mab[i] in b)
  //   match_ba         #E by N matrix, describing for each edge the matching b->a, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #mba[i] of
  //                    the vector set in a is matched to vector #i in b)
  //
  // Output:
  //   singularities    #S by 1 list of the indices of the singular vertices
  //
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedS>
  IGL_INLINE void polyvector_field_singularities_from_matchings(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const std::vector<bool> &V_border,
                                                                const std::vector<std::vector<VFType> > &VF,
                                                                const Eigen::MatrixXi &TT,
                                                                const Eigen::MatrixXi &E2F,
                                                                const Eigen::MatrixXi &F2E,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                Eigen::PlainObjectBase<DerivedS> &singularities);
  
  
  // Wrapper with only V,F and the matchings as input
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedS>
  IGL_INLINE void polyvector_field_singularities_from_matchings(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                Eigen::PlainObjectBase<DerivedS> &singularities);
  
  
  // Same pair as above but also returns singularity indices
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedS>
  IGL_INLINE void polyvector_field_singularities_from_matchings(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                Eigen::PlainObjectBase<DerivedS> &singularity_indices);
  
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedS>
  IGL_INLINE void polyvector_field_singularities_from_matchings(
                                                                const Eigen::PlainObjectBase<DerivedV> &V,
                                                                const Eigen::PlainObjectBase<DerivedF> &F,
                                                                const std::vector<bool> &V_border,
                                                                const std::vector<std::vector<VFType> > &VF,
                                                                const Eigen::MatrixXi &TT,
                                                                const Eigen::MatrixXi &E2F,
                                                                const Eigen::MatrixXi &F2E,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                Eigen::PlainObjectBase<DerivedS> &singularity_indices);

  
  
};


#ifndef IGL_STATIC_LIBRARY
#include "polyvector_field_singularities_from_matchings.cpp"
#endif


#endif /* defined(IGL_POLYVECTOR_FIELD_SINGULARITIES_FROM_MATCHINGS) */
