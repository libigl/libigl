#include "orient_outward_ao.h"
#include "per_face_normals.h"
#include "barycenter.h"
#include "doublearea.h"
#include "matlab_format.h"
#include "embree/ambient_occlusion.h"
#include <iostream>
#include <random>

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedC, 
  typename DerivedFF, 
  typename DerivedI>
IGL_INLINE void igl::orient_outward_ao(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedC> & C,
  const igl::EmbreeIntersector<PointMatrixType,FaceMatrixType,RowVector3> & ei,
  const int num_samples,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace std;
  assert(C.rows() == F.rows());
  assert(F.cols() == 3);
  assert(V.cols() == 3);

  // number of faces
  const int m = F.rows();
  // number of patches
  const int num_cc = C.maxCoeff()+1;
  I.resize(num_cc);
  if(&FF != &F)
  {
    FF = F;
  }
  PlainObjectBase<DerivedV> N;
  Matrix<typename DerivedV::Scalar,Dynamic,1> A;
  per_face_normals(V,F,N);
  doublearea(V,F,A);
  double minarea = A.minCoeff();
  mt19937 engine;
  engine.seed(time(0));
  vector<int> ddist_probability(m);
  for (int f = 0; f < m; ++f)
      ddist_probability[f] = static_cast<int>(A(f) * 100. / minarea);
  discrete_distribution<int> ddist(dist_probability.begin(), dist_probability.end());
  uniform_real_distribution<double> rdist;
  VectorXi face_occluded_front(m, 0);
  VectorXi face_occluded_back (m, 0);
#pragma omp parallel for
  for (int i = 0; i < num_samples; ++i) {
    int f = dist(engine);   // select face with probability proportional to face area
    double t0 = rdist(engine);
    double t1 = rdist(engine);
    double t2 = rdist(engine);
    double t_sum = t0 + t1 + t2;
    t0 /= t_sum;
    t1 /= t_sum;
    t2 /= t_sum;
    RowVector3d p = t0 * V.row(F(f,0)) + t1 * V.row(F(f,1)) + t1 * V.row(F(f,2));
    RowVector3d n = N.row(f);
    bool is_backside = rdist(engine) < 0.5;
    if (is_backside)
        n *= -1;
    Matrix<typename DerivedV::Scalar,Dynamic,1> S;
    ambient_occlusion(ei, p, n, 1, S);
  
  }
  // take area weighted average
  for(int c = 0;c<num_cc;c++)
  {
    //dot(c) /= (typename DerivedV::Scalar) totA(c);
    //if(dot(c) < 0)
    bool b = true;
    I(c) = b;
  }
  // flip according to I
  for(int f = 0;f<m;f++)
  {
    if(I(C(f)))
    {
      FF.row(f) = FF.row(f).reverse().eval();
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::orient_outward_ao<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif

