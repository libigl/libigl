#include "is_border_vertex.h"
#include <vector>

#include "tt.h"

template <typename DerivedV, typename DerivedF>
IGL_INLINE std::vector<bool> igl::is_border_vertex(const Eigen::PlainObjectBase<DerivedV> &V, const Eigen::PlainObjectBase<DerivedF> &F)
{
  Eigen::PlainObjectBase<DerivedF> FF;
  igl::tt(V,F,FF);
  std::vector<bool> ret(V.rows());
  for(unsigned i=0; i<ret.size();++i)
    ret[i] = false;
  
  for(unsigned i=0; i<F.rows();++i)
    for(unsigned j=0;j<F.cols();++j)
      if(FF(i,j) == -1)
      {
        ret[F(i,j)]       = true;
        ret[F(i,(j+1)%F.cols())] = true;
      }
  return ret;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
