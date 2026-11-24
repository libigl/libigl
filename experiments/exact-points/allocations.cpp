#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/point_mesh_closest_point.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>

namespace
{
bool track_allocations=false;
std::size_t allocation_count=0;
}
void *operator new(std::size_t size)
{
  if(track_allocations) { ++allocation_count; }
  if(void *memory=std::malloc(size==0 ? 1 : size)) { return memory; }
  throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory,std::size_t) noexcept { std::free(memory); }

extern "C" int __real_predicates_expansion_sum(int,const double*,int,const double*,double*);
std::size_t sum_calls=0;
extern "C" int __wrap_predicates_expansion_sum(int n,const double *a,int m,const double *b,double *out)
{
  ++sum_calls;
  return __real_predicates_expansion_sum(n,a,m,b,out);
}
int main()
{
  namespace p=igl::predicates;
  using Point=Eigen::RowVector3d;
  const Point a(0,0,0),b(1,0,0),c(0,1,0),q(.25,.25,1);
  int id;
  p::exact_point_triangle_closest_feature(q,a,b,c,id);
  sum_calls=0;
  track_allocations=true;
  const auto dimension=p::exact_point_triangle_closest_feature(q,a,b,c,id);
  const auto filtered_calls=sum_calls;
  const Point d(0,0,1),edge_query(.25,-1,-1);
  const int edge_tie=p::exact_point_triangle_squared_distance_compare(edge_query,a,b,c,b,a,d);
  const auto tie_calls=sum_calls-filtered_calls;
  const double epsilon=std::ldexp(1.0,-52);
  const Point offset(0,0,epsilon),second_b=b+offset,second_c=c+offset,far(0,0,std::ldexp(1.0,50));
  const int order=p::exact_point_triangle_squared_distance_compare(far,a,b,c,offset,second_b,second_c);
  track_allocations=false;
  const auto leaf_allocations=allocation_count,leaf_sums=sum_calls;
  Eigen::MatrixXd vertices(3,3),lower=Eigen::MatrixXd::Zero(1,3),upper=Eigen::MatrixXd::Constant(1,3,2);
  Eigen::MatrixXi faces(1,3); faces<<0,1,2;
  Eigen::VectorXi leaves(1); leaves<<0;
  vertices<<0,0,0, 1,0,0, 0,1,0;
  sum_calls=0; allocation_count=0; track_allocations=true;
  const auto ordinary=p::point_mesh_closest_point(q,vertices,faces,lower,upper,leaves);
  track_allocations=false;
  const auto construction_filtered_calls=sum_calls;
  // This face classification is certified without expansions. Its tiny
  // plane-distance construction must fall back to expansion arithmetic.
  vertices<<0,0,0, 2,2,0, 0,0,2;
  const Point near(.5,.5+epsilon,.5),v0=vertices.row(0),v1=vertices.row(1),v2=vertices.row(2);
  p::exact_point_triangle_closest_feature(near,v0,v1,v2,id);
  const auto construction_feature_calls=sum_calls-construction_filtered_calls;
  track_allocations=true;
  const auto cancelled=p::point_mesh_closest_point(near,vertices,faces,lower,upper,leaves);
  track_allocations=false;
  const auto construction_fallback_calls=sum_calls-construction_filtered_calls-construction_feature_calls;
  std::cout<<"leaf_allocations,filtered_expansion_sums,shared_edge_expansion_sums,fallback_expansion_sums,wrapper_allocations,construction_filtered_sums,construction_feature_sums,construction_fallback_sums\n"
    <<leaf_allocations<<","<<filtered_calls<<","<<tie_calls<<","<<leaf_sums-filtered_calls-tie_calls
    <<","<<allocation_count<<","<<construction_filtered_calls<<","<<construction_feature_calls
    <<","<<construction_fallback_calls<<"\n";
  return dimension!=p::ClosestPointFeature::FACE || order!=1 || edge_tie!=0 ||
    leaf_allocations!=0 || filtered_calls!=0 || tie_calls!=0 || leaf_sums==0 ||
    ordinary.distance!=1 || cancelled.distance<=0 || allocation_count!=2 ||
    construction_filtered_calls!=0 || construction_feature_calls!=0 || construction_fallback_calls==0;
}
