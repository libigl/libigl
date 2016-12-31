#include "redux.h"
#include "for_each.h"

template <typename AType, typename Func, typename DerivedB>
IGL_INLINE void igl::redux(
  const Eigen::SparseMatrix<AType> & A,
  const int dim,
  const Func & func,
  Eigen::PlainObjectBase<DerivedB> & B)
{
  assert((dim == 1 || dim == 2) && "dim must be 2 or 1");
  // Get size of input
  int m = A.rows();
  int n = A.cols();
  // resize output
  B = DerivedB::Zero(dim==1?n:m);
  const auto func_wrap = [&func,&B,&dim](const int i, const int j, const int v)
  {
    if(dim == 1)
    {
      B(j) = i == 0? v : func(B(j),v);
    }else
    {
      B(i) = j == 0? v : func(B(i),v);
    }
  };
  for_each(A,func_wrap);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
