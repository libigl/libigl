#include "vf.h"

#include "verbose.h"

template <typename T, typename S>
IGL_INLINE void igl::vf(
  const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> & V, 
  const Eigen::MatrixXi & F, 
  std::vector<std::vector<T> >& VF, std::vector<std::vector<T> >& VFi)
{
  VF.clear();
  VFi.clear();
  
  VF.resize(V.rows());
  VFi.resize(V.rows());
  
  for(int fi=0; fi<F.rows(); ++fi)
  {
    for(int i = 0; i < 3; ++i)
    {
      VF[F(fi,i)].push_back(fi);
      VFi[F(fi,i)].push_back(i);
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
