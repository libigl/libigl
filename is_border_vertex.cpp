#include "is_border_vertex.h"
#include <vector>

#include "tt.h"

template<typename T, typename S>
IGL_INLINE std::vector<bool> igl::is_border_vertex(const T& V, const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& F)
{
  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> FF;
  igl::tt(V,F,FF);
  std::vector<bool> ret(V.rows());
  for(int i=0; i<ret.size();++i)
    ret[i] = false;
  
  for(int i=0; i<F.rows();++i)
    for(int j=0;j<F.cols();++j)
      if(FF(i,j) == -1)
      {
        ret[F(i,j)]       = true;
        ret[F(i,(j+1)%3)] = true;
      }
  return ret;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
