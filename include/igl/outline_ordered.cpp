#include "outline_ordered.h"

#include "igl/exterior_edges.h"
#include <set>

using namespace std;
using namespace Eigen;

template <typename Index>
IGL_INLINE void igl::outline_ordered(
    const Eigen::MatrixXi& F, 
    std::vector<std::vector<Index> >& L)
{
  MatrixXi E = exterior_edges(F);

  set<int> unseen;
  for (int i = 0; i < E.rows(); ++i)
      unseen.insert(unseen.end(),i);

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
