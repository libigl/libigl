#include "orient_outward_ao.h"
#include "../per_face_normals.h"
#include "../doublearea.h"
#include "../random_dir.h"
#include "EmbreeIntersector.h"
#include <kt84/dbg.h>
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
  
  EmbreeIntersector<float, int> ei(V.cast<float>(),F);
  EmbreeIntersector ei;
  ei.init(V.template cast<float>(),F2.template cast<int>());
  
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
  MatrixXd N;
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
  cout << "generating rays... ";
  uniform_real_distribution<float> rdist;
  mt19937 prng;
  prng.seed(0);
  vector<int     > ray_face;
  vector<Vector3f> ray_ori;
  vector<Vector3f> ray_dir;
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
      float t0 = rdist(prng);        // random barycentric coordinate
      float t1 = rdist(prng);
      float t2 = rdist(prng);
      float t_sum = t0 + t1 + t2;
      t0 /= t_sum;
      t1 /= t_sum;
      t2 /= t_sum;
      Vector3f p = t0 * V.row(F(f,0)).cast<float>().eval()       // be careful with the index!!!
                 + t1 * V.row(F(f,1)).cast<float>().eval()
                 + t2 * V.row(F(f,2)).cast<float>().eval();
      Vector3f n = N.row(f).cast<float>();
      assert(n != Vector3f::Zero());
      // random direction in hemisphere around n (avoid too grazing angle)
      Vector3f d;
      while (true) {
        d = random_dir().cast<float>();
        float ndotd = n.dot(d);
        if (fabsf(ndotd) < 0.1f)
        {
          continue;
        }
        if (ndotd < 0)
        {
          d *= -1.0f;
        }
        break;
      }
      ray_face.push_back(f);
      ray_ori .push_back(p);
      ray_dir .push_back(d);
    }
  }
  
  // per component voting: positive for front side, negative for back side
  vector<float> C_vote_distance(num_cc, 0);     // sum of distance between ray origin and intersection
  vector<int  > C_vote_infinity(num_cc, 0);     // number of rays reaching infinity
  
  auto get_hit_point = [&] (Hit hit) -> Vector3f
  {
    Vector3f p0 = V.row(F(hit.id, 0)).cast<float>();
    Vector3f p1 = V.row(F(hit.id, 1)).cast<float>();
    Vector3f p2 = V.row(F(hit.id, 2)).cast<float>();
    return (1.0f - hit.u - hit.v) * p0 + hit.u * p1 + hit.v * p2;
  };
  
  cout << "shooting rays... ";
#pragma omp parallel for
  for (int i = 0; i < (int)ray_face.size(); ++i)
  {
    int      f = ray_face[i];
    Vector3f o = ray_ori [i];
    Vector3f d = ray_dir [i];
    int c = C(f);
    if (c==65)
      c=c;
    
    // shoot ray toward front & back
    Hit hit_front;
    bool is_hit_front;
    for (float offset = numeric_limits<float>::min(); ; offset *= 10.0f) {
      hit_front.id = -1;
      is_hit_front = ei.intersectRay(o + offset * d, d, hit_front);
      if (!is_hit_front) break;
      if (hit_front.id != f) break;
    }
    Hit hit_back;
    bool is_hit_back;
    for (float offset = numeric_limits<float>::min(); ; offset *= 10.0f) {
      hit_back.id = -1;
      is_hit_back = ei.intersectRay(o - offset * d, -d, hit_back);
      if (!is_hit_back) break;
      if (hit_back.id != f) break;
    }
    if (hit_front.id == f || hit_back.id == f)
    {
      // due to numerical error, ray origin was detected as intersection -> ignore this ray
      continue;
    }
    
    if (is_hit_front)
    {
#pragma omp atomic
      C_vote_distance[c] += (get_hit_point(hit_front) - o).norm();
    } else {
#pragma omp atomic
      C_vote_infinity[c]++;
    }
    
    if (is_hit_back)
    {
#pragma omp atomic
      C_vote_distance[c] -= (get_hit_point(hit_back) - o).norm();
    } else {
#pragma omp atomic
      C_vote_infinity[c]--;
    }
  }
  
  for(int c = 0;c<num_cc;c++)
  {
    I(c) = C_vote_infinity[c] < 0;// || C_vote_distance[c] < 0;
  }
  // flip according to I
  for(int f = 0;f<m;f++)
  {
    if(I(C(f)))
    {
      FF.row(f) = FF.row(f).reverse().eval();
    }
  }
  cout << "done!\n";
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
template void igl::orient_outward_ao<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
