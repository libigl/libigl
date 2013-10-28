#include "ambient_occlusion.h"
#include "EmbreeIntersector.h"
#include <igl/random_dir.h>

template <
  typename PointMatrixType,
  typename FaceMatrixType,
  typename RowVector3,
  typename DerivedP,
  typename DerivedN,
  typename DerivedS >
void igl::ambient_occlusion(
  const igl::EmbreeIntersector<PointMatrixType,FaceMatrixType,RowVector3> & ei,
  const Eigen::PlainObjectBase<DerivedP> & P,
  const Eigen::PlainObjectBase<DerivedN> & N,
  const int num_samples,
  Eigen::PlainObjectBase<DerivedS> & S)
{
  using namespace Eigen;
  using namespace igl;
  const int n = P.rows();
  // Resize output
  S.resize(n,1);
  // Embree seems to be parallel when constructing but not when tracing rays
#pragma omp parallel for
  // loop over mesh vertices
  for(int p = 0;p<n;p++)
  {
    const Vector3d origin = P.row(p);
    const Vector3d normal = N.row(p);
    int num_hits = 0;
    MatrixXd D = random_dir_stratified(num_samples);
    for(int s = 0;s<num_samples;s++)
    {
      //Vector3d d = random_dir();
      Vector3d d = D.row(s);
      if(d.dot(normal) < 0)
      {
        // reverse ray
        d *= -1;
      }
      embree::Hit hit;
      if(ei.intersectRay(origin,d,hit))
      {
        num_hits++;
      }
    }
    S(p) = (double)num_hits/(double)num_samples;
  }
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedP,
  typename DerivedN,
  typename DerivedS >
void igl::ambient_occlusion(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedP> & P,
  const Eigen::PlainObjectBase<DerivedN> & N,
  const int num_samples,
  Eigen::PlainObjectBase<DerivedS> & S)
{
  using namespace igl;
  using namespace Eigen;
  EmbreeIntersector<
    PlainObjectBase<DerivedV>,
    PlainObjectBase<DerivedF>,
    Matrix<typename DerivedP::Scalar,3,1> > ei(V,F);
  return ambient_occlusion(ei,P,N,num_samples,S);
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::ambient_occlusion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::EmbreeIntersector<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, const int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::ambient_occlusion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::ambient_occlusion<Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >, Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(igl::EmbreeIntersector<Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >, Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
