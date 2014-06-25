// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "boundary_loop.h"

#include "triangle_triangle_adjacency.h"
#include "vertex_triangle_adjacency.h"

IGL_INLINE void igl::boundary_loop(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::VectorXi& b)
{
  std::vector<int> bnd;
  bnd.clear();
  std::vector<bool> isVisited(V.rows(),false);

  // Actually mesh only needs to be manifold near boundary, so this is
  // over zealous (see gptoolbox's outline_loop for a more general
  // (and probably faster) implementation)
  assert(is_edge_manifold(V,F) && "Mesh must be manifold");
  Eigen::MatrixXi TT,TTi;
  std::vector<std::vector<int> > VF, VFi;
  igl::triangle_triangle_adjacency(V,F,TT,TTi);
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  // Extract one boundary edge
  bool done = false;
  for (int i = 0; i < TT.rows() && !done; i++)
  {
    for (int j = 0; j < TT.cols(); j++)
    {
      if (TT(i,j) < 0)
      {
        int idx1, idx2;
        idx1 = F(i,j);
        idx2 = F(i,(j+1) % F.cols());
        bnd.push_back(idx1);
        bnd.push_back(idx2);
        isVisited[idx1] = true;
        isVisited[idx2] = true;
        done = true;
        break;
      }
    }
  }

  // Traverse boundary
  while(1)
  {
    bool changed = false;
    int lastV;
    lastV = bnd[bnd.size()-1];

    for (int i = 0; i < (int)VF[lastV].size(); i++)
    {
      int curr_neighbor = VF[lastV][i];

      if (TT.row(curr_neighbor).minCoeff() < 0.) // Face contains boundary edge
      {
        int idx_lastV_in_face;
        if (F(curr_neighbor,0) == lastV) idx_lastV_in_face = 0;
        if (F(curr_neighbor,1) == lastV) idx_lastV_in_face = 1;
        if (F(curr_neighbor,2) == lastV) idx_lastV_in_face = 2;

        int idx_prev = (idx_lastV_in_face + F.cols()-1) % F.cols();
        int idx_next = (idx_lastV_in_face + 1) % F.cols();
        bool isPrevVisited = isVisited[F(curr_neighbor,idx_prev)];
        bool isNextVisited = isVisited[F(curr_neighbor,idx_next)];

        bool gotBndEdge = false;
        int next_bnd_vertex;
        if (!isNextVisited && TT(curr_neighbor,idx_lastV_in_face) < 0)
        {
          next_bnd_vertex = idx_next;
          gotBndEdge = true;
        }
        else if (!isPrevVisited && TT(curr_neighbor,(idx_lastV_in_face+2) % F.cols()) < 0)
        {
          next_bnd_vertex = idx_prev;
          gotBndEdge = true;
        }

        if (gotBndEdge)
        {
          changed = true;
          bnd.push_back(F(curr_neighbor,next_bnd_vertex));
          isVisited[F(curr_neighbor,next_bnd_vertex)] = true;
          break;
        }
      }
    }

    if (!changed)
      break;
  }

  b.resize(bnd.size());
  for(unsigned i=0;i<bnd.size();++i)
    b(i) = bnd[i];
}
