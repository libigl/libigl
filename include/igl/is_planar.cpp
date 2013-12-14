#include "is_planar.h"
IGL_INLINE bool igl::is_planar(const Eigen::MatrixXd & V)
{
  if(V.size() == 0) return false;
  if(V.cols() == 2) return true;
  for(int i = 0;i<V.rows();i++)
  {
    if(V(i,2) != 0) return false;
  }
  return true;
}
