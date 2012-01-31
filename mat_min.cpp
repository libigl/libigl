#include "mat_min.h"

#include "verbose.h"

template <typename T>
IGL_INLINE void igl::mat_min(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
  const int dim,
  Eigen::Matrix<T,Eigen::Dynamic,1> & Y,
  Eigen::Matrix<int,Eigen::Dynamic,1> & I)
{
  assert(dim==1||dim==2);

  // output size
  int n = (dim==1?X.cols():X.rows());
  // resize output
  Y.resize(n);
  I.resize(n);

  // loop over dimension opposite of dim
  for(int j = 0;j<n;j++)
  {
    typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Index PHONY;
    typename Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Index i;
    T m;
    if(dim==1)
    {
      m = X.col(j).minCoeff(&i,&PHONY);
    }else
    {
      m = X.row(j).minCoeff(&PHONY,&i);
    }
    Y(j) = m;
    I(j) = i;
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
