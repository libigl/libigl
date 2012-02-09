#include "matrix_to_list.h"

#include <Eigen/Dense>

template <typename Mat, typename T>
IGL_INLINE void igl::matrix_to_list(
  const Mat & M, 
  std::vector<std::vector<T > > & V)
{
  using namespace std;
  V.resize(M.rows(),vector<T >(M.cols()));
  // loop over rows
  for(int i = 0;i<M.rows();i++)
  {
    // loop over cols
    for(int j = 0;j<M.cols();j++)
    {
      V[i][j] = M(i,j);
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::matrix_to_list<Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >, int>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&);
#endif
