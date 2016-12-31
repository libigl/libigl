#include "for_each.h"

template <typename AType, typename Func>
IGL_INLINE void igl::for_each(
  const Eigen::SparseMatrix<AType> & A,
  const Func & func)
{
  // Can **not** use parallel for because this must be _in order_
  // Iterate over outside
  for(int k=0; k<A.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<AType>::InnerIterator it (A,k); it; ++it)
    {
      func(it.row(),it.col(),it.value());
    }
  }
}

template <typename DerivedA, typename Func>
IGL_INLINE void igl::for_each(
  const Eigen::DenseBase<DerivedA> & A,
  const Func & func)
{
  // Can **not** use parallel for because this must be _in order_
  if(A.IsRowMajor)
  {
    for(typename DerivedA::Index i = 0;i<A.rows();i++)
    {
      for(typename DerivedA::Index j = 0;j<A.cols();j++)
      {
        func(i,j,A(i,j));
      }
    }
  }else
  {
    for(typename DerivedA::Index j = 0;j<A.cols();j++)
    {
      for(typename DerivedA::Index i = 0;i<A.rows();i++)
      {
        func(i,j,A(i,j));
      }
    }
  }
  
}

#include <functional>
