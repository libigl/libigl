#include "exterior_edges.h"
#include "all_edges.h"

#include <cassert>
#include <unordered_map>
#include <utility>

//template <typename T> inline int sgn(T val) {
//      return (T(0) < val) - (val < T(0));
//}

//static void mod2(std::pair<const std::pair<const int, const int>, int>& p)
//{
//  using namespace std;
//  // Be sure that sign of mod matches sign of argument
//  p.second = p.second%2 ? sgn(p.second) : 0;
//}

//// http://stackoverflow.com/a/5517869/148668
//struct Compare
//{
//   int i;
//   Compare(const int& i) : i(i) {}
//};
//bool operator==(const std::pair<std::pair<const int, const int>,int>&p, const Compare& c)
//{
//  return c.i == p.second;
//}
//bool operator==(const Compare& c, const std::pair<std::pair<const int, const int>, int> &p)
//{
//  return c.i == p.second;
//}

IGL_INLINE void igl::exterior_edges(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & E)
{
  using namespace Eigen;
  using namespace std;
  assert(F.cols() == 3);
  MatrixXi all_E;
  all_edges(F,all_E);
  long int n = F.maxCoeff()+1;
  int m = F.minCoeff();
  const auto & compress = [n,m](const int i, const int j)->long int
  {
    return n*(i-m)+(j-m);
  };
  const auto & decompress = [n,m](const long int l,int & i, int & j)
  {
    i = (l / n) + m;
    j = (l % n) + m;
  };

  // Count occurances looking only at pairs (i,j) where i<j, so we count and
  // edge i-->j as +1 and j-->i as -1
  unordered_map<long int,int> C;
  // Loop over edges
  for(int e = 0;e<all_E.rows();e++)
  {
    int i = all_E(e,0);
    int j = all_E(e,1);
    if(i<j)
    {
      // Forward direction --> +1
      C[compress(i,j)]++;
    }else
    {
      // Backward direction --> -1
      C[compress(j,i)]--;
    }
  }
  // Q: Why mod off factors of 2? +1/-1 already takes care of interior edges?
  //// Mod off any factors of 2
  //for_each(C.begin(),C.end(),mod2);
  //int zeros = (int) count(C.begin(),C.end(),Compare(0));
  //E.resize(C.size() - zeros,2);
  E.resize(all_E.rows(),all_E.cols());
  int e = 0;
  // Find all edges with -1 or 1 occurances, flipping direction for -1
  for(const auto & cit : C)
  {
    int i,j;
    if(cit.second > 0)
    {
      decompress(cit.first,i,j);
    } else if(cit.second < 0)
    {
      decompress(cit.first,j,i);
    } else if(cit.second == 0)
    {
      continue;
    }
    for(int k = 0;k<abs(cit.second);k++)
    {
      E(e,0) = i;
      E(e,1) = j;
      e++;
    }
  }
  E.conservativeResize(e,2);
  assert(e == E.rows());
}

IGL_INLINE Eigen::MatrixXi igl::exterior_edges( const Eigen::MatrixXi & F)
{
  using namespace Eigen;
  MatrixXi E;
  exterior_edges(F,E);
  return E;
}
