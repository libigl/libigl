#include "normalize_rows.h"

template <typename DerivedV>
IGL_INLINE void igl::normalize_rows(
                               const Eigen::PlainObjectBase<DerivedV>& A,
                               Eigen::PlainObjectBase<DerivedV> & B)
{
  // Resize output
  B.resize(A.rows(),A.cols());

  // loop over rows
  for(int i = 0; i < A.rows();i++)
  {
    B.row(i) = A.row(i).normalized();
  }
}
