#include "offset_surface.h"
#include "marching_cubes.h"
#include "../voxel_grid.h"
#include <cassert>
#include <iostream>

template <
  typename DerivedV,
  typename DerivedF,
  typename isolevelType,
  typename sType,
  typename DerivedSV,
  typename DerivedSF,
  typename DerivedGV,
  typename Derivedside,
  typename DerivedS>
void igl::copyleft::offset_surface(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const isolevelType isolevel,
  const sType s,
  const SignedDistanceType & signed_distance_type,
  Eigen::PlainObjectBase<DerivedSV> & SV,
  Eigen::PlainObjectBase<DerivedSF> & SF,
  Eigen::PlainObjectBase<DerivedGV> & GV,
  Eigen::PlainObjectBase<Derivedside> & side,
  Eigen::PlainObjectBase<DerivedS> & S)
{
  {
    typedef typename DerivedV::Scalar Scalar;
    Eigen::AlignedBox<Scalar,3> box;
    typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
    assert(V.cols() == 3 && "V must contain positions in 3D");
    RowVector3S min_ext = V.colwise().minCoeff().array() - isolevel;
    RowVector3S max_ext = V.colwise().maxCoeff().array() + isolevel;
    box.extend(min_ext.transpose());
    box.extend(max_ext.transpose());
    igl::voxel_grid(box,s,1,GV,side);
  }
  {
    Eigen::VectorXi I;
    Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,3> C,N;
    igl::signed_distance(
      GV,V,F,signed_distance_type,S,I,C,N);
  }
  DerivedS SS = S.array()-isolevel;
  igl::copyleft::marching_cubes(SS,GV,side(0),side(1),side(2),SV,SF);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::copyleft::offset_surface<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, 1, 3, 1, 1, 3>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, double, int, igl::SignedDistanceType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, 1, 3, 1, 1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
