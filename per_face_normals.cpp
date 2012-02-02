#include "per_face_normals.h"

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::per_face_normals(
                                 const Eigen::PlainObjectBase<DerivedV>& V,
                                 const Eigen::PlainObjectBase<DerivedF>& F,
                                 Eigen::PlainObjectBase<DerivedV> & N)
{
  N.resize(F.rows(),3);
  // loop over faces
  for(int i = 0; i < F.rows();i++)
  {
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = V.row(F(i,1)) - V.row(F(i,0));
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = V.row(F(i,2)) - V.row(F(i,0));
    N.row(i) = (v1.cross(v2)).normalized();
  }
}
