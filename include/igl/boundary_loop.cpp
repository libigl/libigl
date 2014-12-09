// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "boundary_loop.h"
#include "igl/boundary_facets.h"
#include "igl/slice.h"
#include <set>

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
    const Eigen::PlainObjectBase<DerivedF> & F, 
    std::vector<std::vector<Index> >& L)
{
  using namespace std;
  using namespace Eigen;
  MatrixXi E;
  boundary_facets(F, E);

  set<int> unseen;
  for (int i = 0; i < E.rows(); ++i)
  {
    unseen.insert(unseen.end(),i);
  }

  while (!unseen.empty())
  {
      vector<Index> l;

      // Get first vertex of loop
      int startEdge = *unseen.begin();
      unseen.erase(unseen.begin());

      int start = E(startEdge,0);
      int next = E(startEdge,1);
      l.push_back(start);

      while (start != next)
      {
          l.push_back(next);

          // Find next edge
          int nextEdge;
          set<int>::iterator it;
          for (it=unseen.begin(); it != unseen.end() ; ++it)
          {
              if (E(*it,0) == next || E(*it,1) == next)
              {
                  nextEdge = *it;
                  break;
              }                  
          }
          unseen.erase(nextEdge);
          next = (E(nextEdge,0) == next) ? E(nextEdge,1) : E(nextEdge,0);
      }
      L.push_back(l);
  }
}

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
  const Eigen::PlainObjectBase<DerivedF>& F, 
  std::vector<Index>& L)
{
  using namespace Eigen;
  using namespace std;

  vector<vector<int> > Lall;
  boundary_loop(F,Lall);

  int idxMax = -1;
  int maxLen = 0;
  for (int i = 0; i < Lall.size(); ++i)
  {
    if (Lall[i].size() > maxLen)
    {
      maxLen = Lall[i].size();
      idxMax = i;
    }
  }   

  L.resize(Lall[idxMax].size());
  for (int i = 0; i < Lall[idxMax].size(); ++i)
    L[i] = Lall[idxMax][i];
}

template <typename DerivedF, typename DerivedL>
IGL_INLINE void igl::boundary_loop(
  const Eigen::PlainObjectBase<DerivedF>& F, 
  Eigen::PlainObjectBase<DerivedL>& L)
{
  using namespace Eigen;
  using namespace std;

  vector<int> Lvec;
  boundary_loop(F,Lvec);

  L.resize(Lvec.size());
  for (int i = 0; i < Lvec.size(); ++i)
    L(i) = Lvec[i];
}