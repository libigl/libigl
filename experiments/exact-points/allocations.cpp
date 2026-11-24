#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
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
  std::cout<<"leaf_allocations,filtered_expansion_sums,shared_edge_expansion_sums,fallback_expansion_sums\n"
    <<allocation_count<<","<<filtered_calls<<","<<tie_calls<<","<<sum_calls-filtered_calls-tie_calls<<"\n";
  return dimension!=p::ClosestPointFeature::FACE || order!=1 || edge_tie!=0 ||
    allocation_count!=0 || filtered_calls!=0 || tie_calls!=0 || sum_calls==0;
}
