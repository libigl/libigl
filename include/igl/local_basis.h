#ifndef IGL_LOCALBASIS_H
#define IGL_LOCALBASIS_H

#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace igl 
{
  // Compute a local orthogonal reference system for each triangle in the given mesh
  // Templates:
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  // Inputs:
  //   V  eigen matrix #V by 3
  //   F  #F by 3 list of mesh faces (must be triangles)
  // Outputs:
  //   B1 eigen matrix #F by 3, each vector is tangent to the triangle
  //   B2 eigen matrix #F by 3, each vector is tangent to the triangle and perpendicular to B1
  //   B3 eigen matrix #F by 3, normal of the triangle
  //
  // See also: adjacency_matrix
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void local_basis(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV>& B1,
    Eigen::PlainObjectBase<DerivedV>& B2,
    Eigen::PlainObjectBase<DerivedV>& B3
    );

}

#ifdef IGL_HEADER_ONLY
#  include "local_basis.cpp"
#endif

#endif
