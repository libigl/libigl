#include "orient_outward_ao.h"
#include "../per_face_normals.h"
#include "../barycenter.h"
#include "../doublearea.h"
#include "../matlab_format.h"
#include "ambient_occlusion.h"
#include "EmbreeIntersector.h"
#include <iostream>
#include <random>
#include <omp.h>

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedC, 
  typename PointMatrixType,
  typename FaceMatrixType,
  typename RowVector3,
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
  
  // face normal
  PlainObjectBase<DerivedV> N;
  per_face_normals(V,F,N);
  
  // random number generator/distribution for each thread
  int max_threads = omp_get_max_threads();
  
  // prng
  vector<mt19937> engine(max_threads);
  for (int i = 0; i < max_threads; ++i)
      engine[i].seed(time(0) * (i + 1));
  
  // discrete distribution for random selection of faces with probability proportional to their areas
  Matrix<typename DerivedV::Scalar,Dynamic,1> A;
  doublearea(V,F,A);
  double minarea = A.minCoeff();
  Matrix<int, Dynamic, 1> A_int = (A * 100.0 / minarea).template cast<int>();       // only integer is allowed for weight
  auto ddist_func = [&] (double i) { return A_int(static_cast<int>(i)); };
  vector<discrete_distribution<int>> ddist(max_threads, discrete_distribution<int>(m, 0, m, ddist_func));      // simple ctor of (Iter, Iter) not provided by the stupid VC11 impl...
  
  // uniform real between in [0, 1]
  vector<uniform_real_distribution<double>> rdist(max_threads);
  
  // occlusion count per component: +1 when front ray is occluded, -1 when back ray is occluded
  Matrix<int, Dynamic, 1> C_occlude_count;
  C_occlude_count.setZero(m, 1);
  
#pragma omp parallel for
  for (int i = 0; i < num_samples; ++i)
  {
    int thread_num = omp_get_thread_num();
    int f     = ddist[thread_num](engine[thread_num]);   // select face with probability proportional to face area
    double t0 = rdist[thread_num](engine[thread_num]);
    double t1 = rdist[thread_num](engine[thread_num]);
    double t2 = rdist[thread_num](engine[thread_num]);
    double t_sum = t0 + t1 + t2;
    t0 /= t_sum;
    t1 /= t_sum;
    t2 /= t_sum;
    RowVector3d p = t0 * V.row(F(f,0)) + t1 * V.row(F(f,1)) + t1 * V.row(F(f,2));
    RowVector3d n = N.row(f);
    bool is_backside = rdist[thread_num](engine[thread_num]) < 0.5;
    if (is_backside)
    {
        n *= -1;
    }
    Matrix<typename DerivedV::Scalar,Dynamic,1> S;
    ambient_occlusion(ei, p, n, 1, S);
    if (S(0) > 0)
    {
#pragma omp atomic
        C_occlude_count(C(f)) += is_backside ? -1 : 1;
    }

  }
  
  for(int c = 0;c<num_cc;c++)
  {
    I(c) = C_occlude_count(c) > 0;
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

// EmbreeIntersector generated on the fly
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
  const int num_samples,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace igl;
  using namespace Eigen;
  EmbreeIntersector<
    PlainObjectBase<DerivedV>,
    PlainObjectBase<DerivedF>,
    Matrix<typename DerivedV::Scalar,3,1> > ei(V,F);
  return orient_outward_ao(V, F, C, ei, num_samples, FF, I);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::orient_outward_ao<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
