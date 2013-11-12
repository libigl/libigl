#include "moveVF.h"

template <typename T, typename I>
IGL_INLINE void igl::moveVF(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
            const Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic> &F,
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &S,
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &SF)
{
  
  SF = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(F.rows(),S.cols());

  for (int i = 0; i <F.rows(); ++i)
    for (int j = 0; j<F.cols(); ++j)
      SF.row(i) += S.row(F(i,j));

  SF.array() /= F.cols();
  
};

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
