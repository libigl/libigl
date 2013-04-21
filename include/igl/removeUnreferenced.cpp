#include "removeUnreferenced.h"

template <typename T, typename S>
IGL_INLINE void igl::removeUnreferenced(
  const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
  const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> &F,
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV,
  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> &NF,
  Eigen::Matrix<S, Eigen::Dynamic, 1> &I)
{

  // Mark referenced vertices
  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> mark = Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>::Zero(V.rows(),1);
  
  for(int i=0; i<F.rows(); ++i)
  {
    for(int j=0; j<F.cols(); ++j)
    {
      if (F(i,j) != -1)
        mark(F(i,j),0) = 1;
    }
  }
  
  // Sum the occupied cells 
  int newsize = mark.sum();
  
  NV = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(newsize,V.cols());
  NF = Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>(F.rows(),F.cols());
  I  = Eigen::Matrix<S, Eigen::Dynamic, 1>(V.rows(),1);
  
  // Do a pass on the marked vector and remove the unreferenced vertices
  int count = 0;
  for(int i=0;i<mark.rows();++i)
  {
    if (mark(i) == 1)
    {
      NV.row(count) = V.row(i);
      I(i) = count;
      count++;
    }
    else
    {
      I(i) = -1;
    }
  }
  
  // Apply I on F
      
  // Why is this also removing combinatorially degenerate faces?
  count = 0;
  for (int i =0; i<F.rows(); ++i)
  {
    int v0 = I[F(i,0)];
    int v1 = I[F(i,1)];
    int v2 = I[F(i,2)];
    if ( (v0 != v1) && (v2 != v1) && (v0 != v2) )
      NF.row(count++) << v0, v1, v2;
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
