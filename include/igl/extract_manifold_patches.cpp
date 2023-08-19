// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "extract_manifold_patches.h"
#include "unique_edge_map.h"
#include <cassert>
#include <limits>
#include <queue>

template <
  typename DerivedF,
  typename DerivedEMAP,
  typename DeriveduEC,
  typename DeriveduEE,
  typename DerivedP>
IGL_INLINE size_t igl::extract_manifold_patches(
  const Eigen::MatrixBase<DerivedF>& F,
  const Eigen::MatrixBase<DerivedEMAP>& EMAP,
  const Eigen::MatrixBase<DeriveduEC>& uEC,
  const Eigen::MatrixBase<DeriveduEE>& uEE,
  Eigen::PlainObjectBase<DerivedP>& P)
{
  assert(F.cols() == 3);
  const size_t num_faces = F.rows();

  // given directed edge index return corresponding face index
  auto e2f = [&](size_t ei) { return ei % num_faces; };
  // given face index and corner given index of opposite directed edge
  auto fc2e = [&](size_t fi, size_t ci) { return ci*num_faces + fi; };
  // given face index and corner return whether unique edge has incidence == 2
  auto is_manifold_edge = [&](size_t fi, size_t ci) -> bool 
  {
    const size_t ei = fc2e(fi, ci);
    const size_t u = EMAP(ei);
    //return uE2E[EMAP(ei)].size() == 2;
    return (uEC(u+1)-uEC(u))==2;
  };
  // given face and corner index return adjacent face index (assumes manifold)
  auto fc2a = [&](size_t fi, size_t ci) -> size_t 
  {
    const size_t ei = fc2e(fi, ci);
    const size_t u = EMAP(ei);
    //const auto& adj_faces = uE2E[EMAP(ei)];
    //assert(adj_faces.size() == 2);
    assert(uEC(u+1)-uEC(u) == 2);
    //if (adj_faces[0] == ei) 
    //{
    //  return e2f(adj_faces[1]);
    //} else 
    //{
    //  assert(adj_faces[1] == ei);
    //  return e2f(adj_faces[0]);
    //}
    return e2f(uEE(uEC(u)) == ei ? uEE(uEC(u)+1) : uEE(uEC(u)));
  };

  typedef typename DerivedP::Scalar Scalar;
  const Scalar INVALID = std::numeric_limits<Scalar>::max();
  P.resize(num_faces,1);
  P.setConstant(INVALID);
  size_t num_patches = 0;
  for (size_t i=0; i<num_faces; i++) 
  {
    if (P(i,0) != INVALID) continue;

    std::queue<size_t> Q;
    Q.push(i);
    P(i,0) = num_patches;
    while (!Q.empty()) 
    {
      const size_t fid = Q.front();
      Q.pop();
      for (size_t j=0; j<3; j++) 
      {
        if (is_manifold_edge(fid, j)) 
        {
          const size_t adj_fid = fc2a(fid, j);
          if (P(adj_fid,0) == INVALID) 
          {
            Q.push(adj_fid);
            P(adj_fid,0) = num_patches;
          }
        }
      }
    }
    num_patches++;
  }
  assert((P.array() != INVALID).all());

  return num_patches;
}

template<
  typename DerivedF,
  typename DerivedEMAP,
  typename uE2EType,
  typename DerivedP>
IGL_INLINE size_t igl::extract_manifold_patches(
  const Eigen::MatrixBase<DerivedF>& F,
  const Eigen::MatrixBase<DerivedEMAP>& EMAP,
  const std::vector<std::vector<uE2EType> >& uE2E,
  Eigen::PlainObjectBase<DerivedP>& P)
{
    assert(F.cols() == 3);
    const size_t num_faces = F.rows();

    auto edge_index_to_face_index = [&](size_t ei) { return ei % num_faces; };
    auto face_and_corner_index_to_edge_index = [&](size_t fi, size_t ci) {
        return ci*num_faces + fi;
    };
    auto is_manifold_edge = [&](size_t fi, size_t ci) -> bool {
        const size_t ei = face_and_corner_index_to_edge_index(fi, ci);
        return uE2E[EMAP(ei, 0)].size() == 2;
    };
    auto get_adj_face_index = [&](size_t fi, size_t ci) -> size_t {
        const size_t ei = face_and_corner_index_to_edge_index(fi, ci);
        const auto& adj_faces = uE2E[EMAP(ei, 0)];
        assert(adj_faces.size() == 2);
        if (adj_faces[0] == ei) {
            return edge_index_to_face_index(adj_faces[1]);
        } else {
            assert(adj_faces[1] == ei);
            return edge_index_to_face_index(adj_faces[0]);
        }
    };

    typedef typename DerivedP::Scalar Scalar;
    const Scalar INVALID = std::numeric_limits<Scalar>::max();
    P.resize(num_faces,1);
    P.setConstant(INVALID);
    size_t num_patches = 0;
    for (size_t i=0; i<num_faces; i++) {
        if (P(i,0) != INVALID) continue;

        std::queue<size_t> Q;
        Q.push(i);
        P(i,0) = num_patches;
        while (!Q.empty()) {
            const size_t fid = Q.front();
            Q.pop();
            for (size_t j=0; j<3; j++) {
                if (is_manifold_edge(fid, j)) {
                    const size_t adj_fid = get_adj_face_index(fid, j);
                    if (P(adj_fid,0) == INVALID) {
                        Q.push(adj_fid);
                        P(adj_fid,0) = num_patches;
                    }
                }
            }
        }
        num_patches++;
    }
    assert((P.array() != INVALID).all());

    return num_patches;
}

template <
    typename DerivedF,
    typename DerivedP>
IGL_INLINE size_t igl::extract_manifold_patches(
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedP> &P)
{
  Eigen::MatrixXi E, uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<size_t> > uE2E;
  igl::unique_edge_map(F, E, uE, EMAP, uE2E);
  return igl::extract_manifold_patches(F, EMAP, uE2E, P);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template size_t igl::extract_manifold_patches<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template size_t igl::extract_manifold_patches<Eigen::Matrix<int, -1, -1, 0, -1, -1>,   Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template size_t igl::extract_manifold_patches<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, unsigned long, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template size_t igl::extract_manifold_patches<Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, unsigned long, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&); 
#ifdef WIN32
template uint64_t igl::extract_manifold_patches<class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, 1, 0, -1, 1>, uint64_t, class Eigen::Matrix<int, -1, 1, 0, -1, 1>>(class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, class std::vector<class std::vector<uint64_t, class std::allocator<uint64_t>>, class std::allocator<class std::vector<uint64_t, class std::allocator<uint64_t>>>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1>> &);
template uint64_t igl::extract_manifold_patches<class Eigen::Matrix<int, -1, 3, 1, -1, 3>, class Eigen::Matrix<int, -1, 1, 0, -1, 1>, uint64_t, class Eigen::Matrix<int, -1, 1, 0, -1, 1>>(class Eigen::MatrixBase<class Eigen::Matrix<int, -1, 3, 1, -1, 3>> const &, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1>> const &, class std::vector<class std::vector<uint64_t, class std::allocator<uint64_t>>, class std::allocator<class std::vector<uint64_t, class std::allocator<uint64_t>>>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1>> &);
template uint64_t igl::extract_manifold_patches<class Eigen::Matrix<int,-1,-1,0,-1,-1>,class Eigen::Matrix<int,-1,1,0,-1,1>,class Eigen::Matrix<int,-1,1,0,-1,1>,class Eigen::Matrix<int,-1,1,0,-1,1>,class Eigen::Matrix<int,-1,1,0,-1,1> >(class Eigen::MatrixBase<class Eigen::Matrix<int,-1,-1,0,-1,-1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<int,-1,1,0,-1,1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<int,-1,1,0,-1,1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<int,-1,1,0,-1,1> > const &,class Eigen::PlainObjectBase<class Eigen::Matrix<int,-1,1,0,-1,1> > &);
#endif
#endif
