#include <test_common.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/exact_ray_triangle_intersect.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>
#include <exception>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <pthread.h>
#endif

namespace
{
namespace p=igl::predicates;
using Point=Eigen::RowVector3d;
void check_small_stack_queries()
{
  const auto check=[](bool condition)
  { if(!condition) throw std::runtime_error("Small-stack exact query returned an incorrect result"); };
  const double epsilon=std::ldexp(1.,-52);
  const Point a(0,0,0),b(1,0,0),c(0,1,0),q(.25,.25,1);
  const Point offset(0,0,epsilon),second_b=b+offset,second_c=c+offset;
  const Point far(0,0,std::ldexp(1.,50));
  // Cancellation forces exact comparison; coincident faces force an exact tie.
  check(p::exact_point_triangle_squared_distance_compare(far,a,b,c,offset,second_b,second_c)==1);
  check(p::exact_point_triangle_squared_distance_compare(q,a,b,c,c,b,a)==0);
  const Point origin(0,0,0),direction(0,0,1),ra(0,0,1),rb(1,1+epsilon,1),rc(1-epsilon,1,1),
    sa(0,0,2),sb(1,1+epsilon,2),sc(1-epsilon,1,2);
  check(p::exact_ray_triangle_intersect(origin,direction,ra,rb,rc)==p::RayTriangleIntersection::Hit);
  check(p::exact_ray_triangle_hit_compare(origin,direction,ra,rb,rc,sa,sb,sc)==-1);
  check(p::exact_ray_triangle_hit_compare(origin,direction,ra,rb,rc,rc,rb,ra)==0);
  int id;
  check(p::exact_point_triangle_closest_feature(Point(.5,-.5,0),a,Point(1,1,0),Point(0,0,1),id)
    ==p::ClosestPointFeature::VERTEX && id==0);

  Eigen::MatrixXd vertices(4,3),lower=Eigen::MatrixXd::Zero(7,3),upper=Eigen::MatrixXd::Constant(7,3,2);
  vertices<<0,0,0, 2,2,0, 0,0,2, 0,1,0;
  Eigen::MatrixXi faces(4,3),neighbors(4,3),neighbor_edges(4,3);
  faces<<0,1,2, 0,3,1, 0,2,3, 1,3,2;
  Eigen::VectorXi leaves(7); leaves<<-1,-1,-1,0,1,2,3;
  for(int f=0;f<4;++f) for(int i=0;i<3;++i)
    for(int g=0;g<4;++g) for(int j=0;j<3;++j)
      if(faces(f,i)==faces(g,(j+1)%3) && faces(f,(i+1)%3)==faces(g,j))
      { neighbors(f,i)=g; neighbor_edges(f,i)=j; }
  p::MeshFeatureExposure exposure;
  p::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,exposure);
  // The allocation audit verifies that this near-plane construction takes
  // the expansion fallback even though its feature is certified by the filter.
  const Point near(.5,.5+epsilon,.5);
  for(bool cached:{false,true})
  {
    p::PointMeshClosestPoint closest;
    p::PointMeshLocation location;
    const auto evaluate=[&](const Point &query)
    {
      return cached ? p::approximate_signed_distance_exact_sign(query,vertices,faces,
        lower,upper,leaves,exposure,closest,location) :
        p::approximate_signed_distance_exact_sign(query,vertices,faces,
          lower,upper,leaves,neighbors,neighbor_edges,closest,location);
    };
    const double distance=evaluate(near);
    check(location==p::PointMeshLocation::INSIDE && closest.feature==p::ClosestPointFeature::FACE);
    check(std::abs(distance/epsilon+1/std::sqrt(2.))<1e-12);
    check(evaluate(Point(1,1,0))==0 && location==p::PointMeshLocation::ON_SURFACE);
    // Equal edge projections at the closest vertex also refine to expansions.
    check(evaluate(Point(0,-1,-1))>0 && location==p::PointMeshLocation::OUTSIDE);
  }
}
struct StackQuery
{
  std::exception_ptr failure;
  void run()
  {
    try { check_small_stack_queries(); }
    catch(...) { failure=std::current_exception(); }
  }
};
#ifdef _WIN32
DWORD WINAPI stack_query(void *argument)
{ static_cast<StackQuery*>(argument)->run(); return 0; }
#else
void *stack_query(void *argument)
{ static_cast<StackQuery*>(argument)->run(); return nullptr; }
#endif
}

TEST_CASE("predicates::exact_point: 512_KiB_worker_stack","[igl][predicates]")
{
  StackQuery query;
  constexpr std::size_t stack_bytes=512*1024;
#ifdef _WIN32
  const HANDLE thread=CreateThread(nullptr,stack_bytes,stack_query,&query,
    STACK_SIZE_PARAM_IS_A_RESERVATION,nullptr);
  REQUIRE(thread!=nullptr);
  const DWORD waited=WaitForSingleObject(thread,INFINITE);
  const BOOL closed=CloseHandle(thread);
  REQUIRE(waited==WAIT_OBJECT_0);
  REQUIRE(closed!=0);
#else
  pthread_attr_t attributes;
  REQUIRE(pthread_attr_init(&attributes)==0);
  const int sized=pthread_attr_setstacksize(&attributes,stack_bytes);
  pthread_t thread;
  const int created=sized==0 ? pthread_create(&thread,&attributes,stack_query,&query) : sized;
  const int destroyed=pthread_attr_destroy(&attributes);
  REQUIRE(created==0);
  const int joined=pthread_join(thread,nullptr);
  REQUIRE(destroyed==0);
  REQUIRE(joined==0);
#endif
  if(query.failure) std::rethrow_exception(query.failure);
}
