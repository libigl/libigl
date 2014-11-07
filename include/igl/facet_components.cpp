#include "facet_components.h"
#include <igl/triangle_triangle_adjacency.h>
#include <vector>
#include <queue>
template <typename DerivedF, typename DerivedC>
IGL_INLINE void igl::facet_components(
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  using namespace std;
  using namespace igl;
  typedef typename DerivedF::Index Index;
  vector<vector<vector<Index > > > TT;
  vector<vector<vector<Index > > > TTi;
  triangle_triangle_adjacency(F,TT,TTi);
  return facet_components(TT,C);
}

template <
  typename TTIndex, 
  typename DerivedC>
IGL_INLINE void igl::facet_components(
  const std::vector<std::vector<std::vector<TTIndex > > > & TT,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  using namespace std;
  using namespace igl;
  typedef TTIndex Index;
  const Index m = TT.size();
  C.resize(m,1);
  vector<bool> seen(m,false);
  Index id = 0;
  for(Index g = 0;g<m;g++)
  {
    if(seen[g])
    {
      continue;
    }
    queue<Index> Q;
    Q.push(g);
    while(!Q.empty())
    {
      const Index f = Q.front();
      Q.pop();
      if(seen[f])
      {
        continue;
      }
      seen[f] = true;
      C(f,0) = id;
      // Face f's neighbor lists opposite opposite each corner
      for(const auto & c : TT[f])
      {
        // Each neighbor
        for(const auto & n : c)
        {
          if(!seen[n])
          {
            Q.push(n);
          }
        }
      }
    }
    id++;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::facet_components<long, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(std::vector<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > >, std::allocator<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
#endif
