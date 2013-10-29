#include "orient_outward_ao.h"
#include "../per_face_normals.h"
#include "../doublearea.h"
#include "../random_dir.h"
#include "EmbreeIntersector.h"
#include <iostream>
#include <random>
#include <limits>

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
  const int min_num_rays_per_component,
  const int total_num_rays,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace std;
  assert(C.rows() == F.rows());
  assert(F.cols() == 3);
  assert(V.cols() == 3);
  
  // pass both sides of faces to Embree
  MatrixXi F2;
  F2.resize(F.rows()*2,F.cols());
  F2 << F, F.rowwise().reverse().eval();
  EmbreeIntersector<typename DerivedV::Scalar, typename DerivedF::Scalar> ei(V,F2);
  
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
  
  // face area
  Matrix<typename DerivedV::Scalar,Dynamic,1> A;
  doublearea(V,F,A);
  double area_min = A.minCoeff();
  double area_total = A.sum();
  
  // determine number of rays per component according to its area
  VectorXd area_per_component;
  area_per_component.setZero(num_cc);
  for (int f = 0; f < m; ++f)
  {
    area_per_component(C(f)) += A(f);
  }
  VectorXi num_rays_per_component;
  num_rays_per_component.setZero(num_cc);
  for (int c = 0; c < num_cc; ++c)
  {
    num_rays_per_component(c) = max<int>(min_num_rays_per_component, static_cast<int>(total_num_rays * area_per_component(c) / area_total));
  }
  
  // generate all the rays
  uniform_real_distribution<double> rdist;
  mt19937 prng;
  prng.seed(time(0));
  vector<int     > ray_face;
  vector<Vector3d> ray_ori;
  vector<Vector3d> ray_dir;
  ray_face.reserve(total_num_rays);
  ray_ori .reserve(total_num_rays);
  ray_dir .reserve(total_num_rays);
  for (int c = 0; c < num_cc; ++c)
  {
    vector<int> CF;     // set of faces per component
    vector<int> CF_area;
    for (int f = 0; f < m; ++f)
    {
      if (C(f)==c)
      {
        CF.push_back(f);
        CF_area.push_back(static_cast<int>(100 * A(f) / area_min));
      }
    }
    // discrete distribution for random selection of faces with probability proportional to their areas
    auto ddist_func = [&] (double i) { return CF_area[static_cast<int>(i)]; };
    discrete_distribution<int> ddist(CF.size(), 0, CF.size(), ddist_func);      // simple ctor of (Iter, Iter) not provided by the stupid VC11 impl...
    for (int i = 0; i < num_rays_per_component[c]; ++i)
    {
      int f     = CF[ddist(prng)];    // select face with probability proportional to face area
      double t0 = rdist(prng);        // random barycentric coordinate
      double t1 = rdist(prng);
      double t2 = rdist(prng);
      double t_sum = t0 + t1 + t2;
      t0 /= t_sum;
      t1 /= t_sum;
      t2 /= t_sum;
      Vector3d p = t0 * V.row(F(f,0))       // be careful with the index!!!
                 + t1 * V.row(F(f,1))
                 + t2 * V.row(F(f,2));
      Vector3d n = N.row(f);
      Vector3d d = random_dir();
      if (n.dot(d) < 0)
      {
        d *= -1;
      }
      ray_face.push_back(f);
      ray_ori .push_back(p);
      ray_dir .push_back(d);
    }
  }
  
  // occlusion count per component
  vector<int> C_occlude_count_front(num_cc, 0);
  vector<int> C_occlude_count_back (num_cc, 0);

  //auto dbg_get_hit_point = [&] (embree::Hit hit) {
  //  RowVector3d p0 = V.row(F2(hit.id0, 0));
  //  RowVector3d p1 = V.row(F2(hit.id0, 1));
  //  RowVector3d p2 = V.row(F2(hit.id0, 2));
  //  RowVector3d p = (1 - hit.u - hit.v) * p0 + hit.u * p1 + hit.v * p2;
  //  return p;
  //};

#pragma omp parallel for
  for (int i = 0; i < ray_face.size(); ++i)
  {
    int      f = ray_face[i];
    Vector3d o = ray_ori [i];
    Vector3d d = ray_dir [i];
    int c = C(f);
    Hit hit_front;
    if (ei.intersectRay(o, d, hit_front))
    {
#pragma omp atomic
      //RowVector3d hit_point = dbg_get_hit_point(hit_front);
      C_occlude_count_front[c]++;
    }
    Hit hit_back;
    if (ei.intersectRay(o, -d, hit_back))
    {
#pragma omp atomic
      //RowVector3d hit_point = dbg_get_hit_point(hit_back);
      C_occlude_count_back[c]++;
    }
  }
  
  for(int c = 0;c<num_cc;c++)
  {
    I(c) = C_occlude_count_front[c] > C_occlude_count_back[c];
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

// Call with default parameters
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
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  return orient_outward_ao(V, F, C, 100, F.rows() * 100, FF, I);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::orient_outward_ao<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
