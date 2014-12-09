// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "outline_ordered.h"
#include "exterior_edges.h"
#include <set>

template <typename DerivedF, typename Index>
IGL_INLINE void igl::outline_ordered(
    const Eigen::PlainObjectBase<DerivedF> & F, 
    std::vector<std::vector<Index> >& L)
{
  using namespace std;
  using namespace Eigen;
  // Exterior edges include some non-manifold edges (poor function name). I
  // suppose `outline_ordered` is not well defined for non-manifold meshes, but
  // perhaps this should just call `boundary_facets`
  MatrixXi E = exterior_edges(F);

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
