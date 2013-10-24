#ifndef IGL_UNIQUE_H
#define IGL_UNIQUE_H
#include "igl_inline.h"

#include <vector>
#include <Eigen/Core>
namespace igl
{
  // Act like matlab's [C,IA,IC] = unique(X)
  //
  // Templates:
  //   T  comparable type T
  // Inputs:
  //   A  #A vector of type T
  // Outputs:
  //   C  #C vector of unique entries in A
  //   IA  #C index vector so that C = A(IA);
  //   IC  #A index vector so that A = C(IC);
  template <typename T>
  IGL_INLINE void unique(
    const std::vector<T> & A,
    std::vector<T> & C,
    std::vector<size_t> & IA,
    std::vector<size_t> & IC);

  // Act like matlab's [C,IA,IC] = unique(X,'rows')
  //
  // Templates:
  //   DerivedA derived scalar type, e.g. MatrixXi or MatrixXd
  //   DerivedIA derived integer type, e.g. MatrixXi
  //   DerivedIC derived integer type, e.g. MatrixXi
  // Inputs:
  //   A  m by n matrix whose entries are to unique'd according to rows
  // Outputs:
  //   C  #C vector of unique rows in A
  //   IA  #C index vector so that C = A(IA,:);
  //   IC  #A index vector so that A = C(IC,:);
  template <typename DerivedA, typename DerivedIA, typename DerivedIC>
  IGL_INLINE void unique_rows(
    const Eigen::PlainObjectBase<DerivedA>& A,
    Eigen::PlainObjectBase<DerivedA>& C,
    Eigen::PlainObjectBase<DerivedIA>& IA,
    Eigen::PlainObjectBase<DerivedIC>& IC);
}

#ifdef IGL_HEADER_ONLY
#  include "unique.cpp"
#endif

#endif


