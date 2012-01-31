#include "is_manifold.h"

#include <algorithm>

template<typename T>
IGL_INLINE bool igl::is_manifold(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F)
  {
    std::vector<std::vector<int> > TTT;
    for(int f=0;f<F.rows();++f)
      for (int i=0;i<3;++i)
      {
        // v1 v2 f ei 
        int v1 = F(f,i);
        int v2 = F(f,(i+1)%3);
        if (v1 > v2) std::swap(v1,v2);
        std::vector<int> r(4);
        r[0] = v1; r[1] = v2;
        r[2] = f;  r[3] = i;
        TTT.push_back(r);
      }
    std::sort(TTT.begin(),TTT.end());
    
    for(int i=2;i<(int)TTT.size();++i)
    {
      std::vector<int>& r1 = TTT[i-2];
      std::vector<int>& r2 = TTT[i-1];
      std::vector<int>& r3 = TTT[i];
      if ( (r1[0] == r2[0] && r2[0] == r3[0]) 
          && 
          (r1[1] == r2[1] && r2[1] == r3[1]) )
      {
        return false;
      }
    }
    return true;
  }

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
