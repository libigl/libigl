#ifndef IGL_BOUNDARY_FACES_H
#define IGL_BOUNDARY_FACES_H
#include "igl_inline.h"

#ifndef IGL_NO_EIGEN
#  include <Eigen/Dense>
#endif

#include <vector>

namespace igl
{
  // BOUNDARY_FACES Determine boundary faces of tetrahedra stored in T
  //
  // Templates:
  //   IntegerT  integer-value: i.e. int
  //   IntegerF  integer-value: i.e. int
  // Input:
  //  T  tetrahedron index list, m by 4, where m is the number of tetrahedra
  // Output:
  //  F  list of boundary faces, n by 3, where n is the number of boundary faces
  //
  // Input:
  //  V  # vertices by 3 vertex positions
  //  F  # faces by 3 list of face indices
  // Output: 
  //  RV  # vertices by 3 vertex positions, order such that if the jth vertex is
  //    some face in F, and the kth vertex is not then j comes before k
  //  RF  # faces by 3 list of face indices, reindexed to use RV
  //  IM  # faces by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  //    and RV(IM,:) = V
  //
  template <typename IntegerT, typename IntegerF>
  IGL_INLINE void boundary_faces(
    const std::vector<std::vector<IntegerT> > & T,
    std::vector<std::vector<IntegerF> > & F);

#ifndef IGL_NO_EIGEN
  // Templates:
  //   DerivedT  integer-value: i.e. from MatrixXi
  //   DerivedF  integer-value: i.e. from MatrixXi
  template <typename DerivedT, typename DerivedF>
  IGL_INLINE void boundary_faces(
    const Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F);
#endif
}

#ifdef IGL_HEADER_ONLY
#  include "boundary_faces.cpp"
#endif

#endif

