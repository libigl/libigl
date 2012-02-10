#include "limit_faces.h"

#include <vector>
#include <Eigen/Dense>

template <typename MatF, typename VecL>
IGL_INLINE void igl::limit_faces(
  const MatF & F, 
  const VecL & L, 
  const bool exclusive,
  MatF & LF)
{
  using namespace std;
  using namespace Eigen;
  vector<bool> in(F.rows(),false);
  int num_in = 0;
  // loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    bool all = true;
    bool any = false;
    for(int j = 0;j<F.cols();j++)
    {
      bool found = false;
      // loop over L
      for(int l = 0;l<L.size();l++)
      {
        if(F(i,j) == L(l))
        {
          found = true;
          break;
        }
      }
      any |= found;
      all &= found;
    }
    in[i] = (exclusive?all:any);
    num_in += (in[i]?1:0);
  }

  LF.resize(num_in,F.cols());
  // loop over faces
  int lfi = 0;
  for(int i = 0;i<F.rows();i++)
  {
    if(in[i])
    {
      LF.row(lfi) = F.row(i);
      lfi++;
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
