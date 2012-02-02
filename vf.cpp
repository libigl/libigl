#include "vf.h"

#include "verbose.h"

template <typename DerivedV, typename DerivedF, typename IndexType>
IGL_INLINE void igl::vf(
                   const Eigen::PlainObjectBase<DerivedV>& V,
                   const Eigen::PlainObjectBase<DerivedF>& F,
                   std::vector<std::vector<IndexType> >& VF,
                   std::vector<std::vector<IndexType> >& VFi)
{
  VF.clear();
  VFi.clear();
  
  VF.resize(V.rows());
  VFi.resize(V.rows());
  
  for(int fi=0; fi<F.rows(); ++fi)
  {
    for(int i = 0; i < F.cols(); ++i)
    {
      VF[F(fi,i)].push_back(fi);
      VFi[F(fi,i)].push_back(i);
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
