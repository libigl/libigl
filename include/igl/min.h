#ifndef IGL_MIN_H
#define IGL_MIN_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  /// @param[in] X  m by n matrix
  /// @param[in] dim  dimension along which to take min
  /// @param[out] Y
  ///   Y  n-long vector (if dim == 1) 
  ///   or
  ///   Y  m-long vector (if dim == 2)
  ///   I  vector the same size as Y containing the indices along dim of minimum
  ///     entries
  ///  \deprecated seems like a duplicate of mat_min
  template <typename AType, typename DerivedB, typename DerivedI>
  IGL_INLINE void min(
    const Eigen::SparseMatrix<AType> & A,
    const int dim,
    Eigen::PlainObjectBase<DerivedB> & B,
    Eigen::PlainObjectBase<DerivedI> & I);
}
#ifndef IGL_STATIC_LIBRARY
#  include "min.cpp"
#endif
#endif

