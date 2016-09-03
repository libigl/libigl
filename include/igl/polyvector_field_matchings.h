// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_POLYVECTOR_FIELD_MATCHINGS
#define IGL_POLYVECTOR_FIELD_MATCHINGS
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl {
  // Given a pair of adjacent faces a and b, and a set of N unordered vectors
  // of a vector set defined on each face, the function computes the best order-preserving
  // matching between them, where "best" means "minimum curl", or "smoothest".
  // Inputs:
  //   _ua              1 by 3N row vector of the vectors on face a, stacked horizontally
  //   _ub              1 by 3N row vector of the vectors on face b, stacked horizontally
  //   na               1 by 3, normal of face a
  //   nb               1 by 3, normal of face b
  //   e                1 by 3, the vector corresponding to the shared edge between a and b
  //   match_with_curl  boolean flag, determines whether a curl or a smoothness matching will
  //                    be computed
  //   is_symmetric     boolean flag, determines whether the input vector set field is symmetric(
  //                    =consists of pairwise collinear vectors in each set, in which case only one
  //                    of the vectors in the pair is stored) or not, i.e. the set contains all the vectors
  // )
  // Outputs:
  //   mab              1 by N row vector, describing the matching a->b (i.e. vector #i of the
  //                    vector set in a is matched to vector #mab[i] in b)
  //   mba              1 by N row vector, describing the matching b->a (i.e. vector #mba[i]
  //                    of the vector set in a is matched to vector #i in b)
  //
  template <typename DerivedS, typename DerivedV, typename DerivedM>
  IGL_INLINE void polyvector_field_matching(
                                            const Eigen::PlainObjectBase<DerivedS>& _ua,
                                            const Eigen::PlainObjectBase<DerivedS>& _ub,
                                            const Eigen::PlainObjectBase<DerivedV>& na,
                                            const Eigen::PlainObjectBase<DerivedV>& nb,
                                            const Eigen::PlainObjectBase<DerivedV>& e,
                                            bool match_with_curl,
                                            Eigen::PlainObjectBase<DerivedM>& mab,
                                            Eigen::PlainObjectBase<DerivedM>& mba,
                                            bool is_symmetric);


  // Given a mesh and a vector set field consisting of unordered N-vector sets defined
  // on the faces of the mesh, the function computes, for each (non-boundary) edge
  // the best order-preserving matching between the vector sets of the faces across
  // the edge, where "best" means to "with minimum curl", or "smoothest"
  // Inputs:
  //   sol3D            #F by 3n list of the 3D coordinates of the per-face vectors
  //                    (stacked horizontally for each triangle)
  //   V                #V by 3 list of mesh vertex coordinates
  //   F                #F by 3 list of mesh faces
  //   E                #E by 2 list of mesh edges (pairs of vertex indices)
  //   FN               #F by 3 list of face normals
  //   E2F              #E by 2 list of the edge-to-face relation (e.g. computed
  //                    via igl::edge_topology)
  //   match_with_curl  boolean flag, determines whether curl or smoothness matchings will
  //                    be computed
  //   is_symmetric     boolean flag, determines whether the input vector set field is symmetric(
  //                    =consists of pairwise collinear vectors in each set, in which case only one
  //                    of the vectors in the pair is stored) or not, i.e. the set contains all the vectors
  // Outputs:
  //   match_ab         #E by N matrix, describing for each edge the matching a->b, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #i of
  //                    the vector set in a is matched to vector #mab[i] in b)
  //   match_ba         #E by N matrix, describing for each edge the matching b->a, where a
  //                    and b are the faces adjacent to the edge (i.e. vector #mba[i] of
  //                    the vector set in a is matched to vector #i in b)
  //   curl             the per-edge curl of the matchings (zero for boundary edges)
  // Returns:
  //   meanCurl         the average of the per-edge curl values (only non-boundary edges are counted)
  //
  template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedM, typename DerivedC>
  IGL_INLINE typename DerivedC::Scalar polyvector_field_matchings(
                                                                  const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                                  const Eigen::PlainObjectBase<DerivedV>&V,
                                                                  const Eigen::PlainObjectBase<DerivedF>&F,
                                                                  const Eigen::PlainObjectBase<DerivedE>&E,
                                                                  const Eigen::PlainObjectBase<DerivedV>& FN,
                                                                  const Eigen::MatrixXi &E2F,
                                                                  bool match_with_curl,
                                                                  bool is_symmetric,
                                                                  Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                                  Eigen::PlainObjectBase<DerivedM>& match_ba,
                                                                  Eigen::PlainObjectBase<DerivedC>& curl);


  //Wrapper of the above with only vertices and faces as mesh input
  template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedC>
  IGL_INLINE typename DerivedC::Scalar polyvector_field_matchings(
                                                                  const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                                  const Eigen::PlainObjectBase<DerivedV>&V,
                                                                  const Eigen::PlainObjectBase<DerivedF>&F,
                                                                  bool match_with_curl,
                                                                  bool is_symmetric,
                                                                  Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                                  Eigen::PlainObjectBase<DerivedM>& match_ba,
                                                                  Eigen::PlainObjectBase<DerivedC>& curl);

  //Wrappers with no curl output
  template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedM>
  IGL_INLINE void polyvector_field_matchings(
                                             const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                             const Eigen::PlainObjectBase<DerivedV>&V,
                                             const Eigen::PlainObjectBase<DerivedF>&F,
                                             bool match_with_curl,
                                             bool is_symmetric,
                                             Eigen::PlainObjectBase<DerivedM>& match_ab,
                                             Eigen::PlainObjectBase<DerivedM>& match_ba);
  template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedM>
  IGL_INLINE void polyvector_field_matchings(
                                             const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                             const Eigen::PlainObjectBase<DerivedV>&V,
                                             const Eigen::PlainObjectBase<DerivedF>&F,
                                             const Eigen::PlainObjectBase<DerivedE>&E,
                                             const Eigen::PlainObjectBase<DerivedV>& FN,
                                             const Eigen::MatrixXi &E2F,
                                             bool match_with_curl,
                                             bool is_symmetric,
                                             Eigen::PlainObjectBase<DerivedM>& match_ab,
                                             Eigen::PlainObjectBase<DerivedM>& match_ba);
  
};


#ifndef IGL_STATIC_LIBRARY
#include "polyvector_field_matchings.cpp"
#endif


#endif /* defined(IGL_POLYVECTOR_FIELD_MATCHINGS) */
