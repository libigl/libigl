// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <zhongshi@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "tetrahedron_tetrahedron_adjacency.h"
#include "all_edges.h"
#include "unique_simplices.h"
#include "parallel_for.h"
#include "unique_edge_map.h"
#include "slice.h"
#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <array>
namespace igl
{
  struct tetrahedron_tetrahedron_adjacency_tet_cell
  {
    int v1;
    int v2;
    int v3;
    int t;
    int f;
  };

// Extract the face adjacencies
  template<typename DerivedT, typename TTT_type, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tetrahedron_tetrahedron_adjacency_extractTTi(
      const Eigen::MatrixBase<DerivedT> &F,
      std::vector<TTT_type> &TTT,
      Eigen::MatrixBase<DerivedTT> &TT,
      Eigen::MatrixBase<DerivedTTif> &TTif,
      Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    TT.derived().setConstant((int) (F.rows()), 4, -1);
    TTif.derived().setConstant((int) (F.rows()), 4, -1);
    TTie.derived().setConstant((int) (F.rows()), 4 * 3, -1);

    for(int i = 1; i < (int) TTT.size(); ++i)
    {
      auto &r1 = TTT[i - 1];
      auto &r2 = TTT[i];

      // if get same face.
      if((r1.v1 == r2.v1) && (r1.v2 == r2.v2) && (r1.v3 == r2.v3))
      {
        TT(r1.t, r1.f) = r2.t;
        TT(r2.t, r2.f) = r1.t;

        // inverse for face
        TTif(r1.t, r1.f) = r2.f;
        TTif(r2.t, r2.f) = r1.f;

        // inverse for edge
        Eigen::RowVector3i f1v, f2v;
        for(int i = 0; i < 3; i++)
        {
          f1v[i] = F(r1.t, tetrahedron_local_FF(r1.f, i));
          f2v[i] = F(r2.t, tetrahedron_local_FF(r2.f, i));
        }

        std::vector<int> emap = {0, 1, 2};
        for(int e = 0; e < 3; e++)
          for(int i = 0; i < 3 - e; i++)
            if(f1v(e) == f2v(emap[i]))
            {
              TTie(r1.t, 3 * r1.f + (e + 2) % 3) = emap[i];
              TTie(r2.t, 3 * r2.f + (emap[i] + 2) % 3) = e;
              emap[i] = emap[2 - e]; // fill with the last
              break;
            }

      } //if
    }
  }

  template<typename DerivedT, typename TTT_type>
  IGL_INLINE void tetrahedron_tetrahedron_adjacency_preprocess(
      const Eigen::MatrixBase<DerivedT> &F,
      std::vector<TTT_type> &TTT)
  {
    for(int t = 0; t < (int) F.rows(); ++t)
      for(int f = 0; f < 4; ++f)
      {
        std::array<int, 3> v;
        for(int i = 0; i < 3; ++i)
        {
          v[i] = F(t, tetrahedron_local_FF(f, i));
        }
        std::sort(v.begin(), v.end());
        TTT.push_back({v[0], v[1], v[2], t, f});
      }
    std::sort(TTT.begin(), TTT.end(),
              [](const tetrahedron_tetrahedron_adjacency_tet_cell &r1,
                 const tetrahedron_tetrahedron_adjacency_tet_cell &r2)
              {
                return std::tie(r1.v1, r1.v2, r1.v3, r1.t) <
                       std::tie(r2.v1, r2.v2, r2.v3, r2.t);
              });
  }
}

// Compute tetrahedron-tetrahedron adjacency with indices
template<typename DerivedT, typename DerivedTT,
         typename DerivedTTif,typename DerivedTTie>
IGL_INLINE void
igl::tetrahedron_tetrahedron_adjacency(const Eigen::MatrixBase<DerivedT> &F,
                                       Eigen::MatrixBase<DerivedTT> &TT,
                                       Eigen::MatrixBase<DerivedTTif> &TTif,
                                       Eigen::MatrixBase<DerivedTTie> &TTie)
{
  std::vector<igl::tetrahedron_tetrahedron_adjacency_tet_cell> TTT;
  tetrahedron_tetrahedron_adjacency_preprocess(F, TTT);
  tetrahedron_tetrahedron_adjacency_extractTTi(F, TTT, TT, TTif, TTie);
}

template <typename DerivedT, typename DerivedTT>
IGL_INLINE void igl::tetrahedron_tetrahedron_adjacency(
    const Eigen::MatrixBase<DerivedT>& F,
    Eigen::MatrixBase<DerivedTT>& TT)
{
  std::vector<igl::tetrahedron_tetrahedron_adjacency_tet_cell> TTT;
  tetrahedron_tetrahedron_adjacency_preprocess(F, TTT);
  Eigen::MatrixBase<DerivedTT> TTif,TTie;
  tetrahedron_tetrahedron_adjacency_extractTTi(F,TTT,TT,TTif,TTie);
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::tetrahedron_tetrahedron_adjacency<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
