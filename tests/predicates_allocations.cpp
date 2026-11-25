#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/point_mesh_closest_point.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/triangle_triangle_adjacency.h>
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
extern "C" double __real_nextafter(double,double);
std::size_t rounding_calls=0;
extern "C" double __wrap_nextafter(double value,double direction)
{
  ++rounding_calls;
  return __real_nextafter(value,direction);
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
  const auto wrapper_allocations=allocation_count;
  // Equal boxes visit the right leaf first. Replacing face 1 by face 0 on
  // an exact tie must reuse the bound, without any further interval rounding.
  vertices<<0,0,0, 1,0,0, 0,1,0;
  const Point vertex_query(-1,-2,3);
  rounding_calls=0;
  const auto single=p::exact_point_mesh_closest_triangle(vertex_query,vertices,faces,lower,upper,leaves);
  const auto single_rounding=rounding_calls;
  Eigen::MatrixXi duplicate_faces(2,3); duplicate_faces<<0,1,2, 0,1,2;
  Eigen::MatrixXd duplicate_lower=lower.replicate(3,1),duplicate_upper=upper.replicate(3,1);
  Eigen::VectorXi duplicate_leaves(3); duplicate_leaves<<-1,0,1;
  rounding_calls=0;
  const auto duplicate=p::exact_point_mesh_closest_triangle(vertex_query,vertices,duplicate_faces,
    duplicate_lower,duplicate_upper,duplicate_leaves);
  const auto tie_rounding=rounding_calls;
  // Distinct edge projections need no expansions. Comparing the seed edge
  // with itself would introduce a redundant exact-zero fallback.
  Eigen::MatrixXd tetra_vertices(4,3);
  tetra_vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
  Eigen::MatrixXi tetra_faces(4,3),neighbors,neighbor_edges;
  tetra_faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
  igl::triangle_triangle_adjacency(tetra_faces,neighbors,neighbor_edges);
  sum_calls=0;
  const auto location=p::exact_point_mesh_classify(Point(-1,-2,-3),tetra_vertices,
    tetra_faces,0,neighbors,neighbor_edges);
  const auto vertex_sign_calls=sum_calls;
  const Point ray_origin(0,0,0),ray_direction(0,0,1),
    first_a(-1,-1,1),first_b(1,-1,1),first_c(0,1,1),
    last_a(-1,-1,3),last_b(1,-1,3),last_c(0,1,3);
  sum_calls=0; allocation_count=0; track_allocations=true;
  const int ray_order=p::exact_ray_triangle_hit_compare(ray_origin,ray_direction,
    first_a,first_b,first_c,last_a,last_b,last_c);
  const auto ray_filtered_calls=sum_calls;
  const int ray_tie=p::exact_ray_triangle_hit_compare(ray_origin,ray_direction,
    first_a,first_b,first_c,first_c,first_b,first_a);
  track_allocations=false;
  const auto ray_allocations=allocation_count,ray_fallback_calls=sum_calls-ray_filtered_calls;
  std::cout<<"leaf_allocations,filtered_expansion_sums,shared_edge_expansion_sums,fallback_expansion_sums,wrapper_allocations,construction_filtered_sums,construction_feature_sums,construction_fallback_sums\n"
    <<leaf_allocations<<","<<filtered_calls<<","<<tie_calls<<","<<leaf_sums-filtered_calls-tie_calls
    <<","<<wrapper_allocations<<","<<construction_filtered_calls<<","<<construction_feature_calls
    <<","<<construction_fallback_calls<<"\n";
  std::cout<<"single_bound_rounding,tied_bound_rounding,vertex_sign_expansion_sums\n"
    <<single_rounding<<","<<tie_rounding<<","<<vertex_sign_calls<<"\n";
  return dimension!=p::ClosestPointFeature::FACE || order!=1 || edge_tie!=0 ||
    leaf_allocations!=0 || filtered_calls!=0 || tie_calls!=0 || leaf_sums==0 ||
    ordinary.distance!=1 || cancelled.distance<=0 ||
    construction_filtered_calls!=0 || construction_feature_calls!=0 || construction_fallback_calls==0 ||
    single!=0 || duplicate!=0 || single_rounding!=3 || tie_rounding!=single_rounding ||
    location!=p::PointMeshLocation::OUTSIDE || vertex_sign_calls!=0 ||
    ray_order!=-1 || ray_tie!=0 || ray_filtered_calls!=0 || ray_fallback_calls==0 || ray_allocations!=0;
}
