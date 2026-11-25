#include <test_common.h>
#include <igl/predicates/exact_ray_triangle_intersect.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>
#include <igl/predicates/exactinit.h>
#include <predicates.h>
#include <array>
#include <cstdint>
#include <random>

namespace
{
using IntegerPoint = std::array<std::int64_t,3>;
IntegerPoint subtract(const IntegerPoint &a, const IntegerPoint &b)
{
  return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}};
}
std::int64_t determinant(const IntegerPoint &a, const IntegerPoint &b, const IntegerPoint &c)
{
  return a[0]*(b[1]*c[2]-b[2]*c[1])
       + a[1]*(b[2]*c[0]-b[0]*c[2])
       + a[2]*(b[0]*c[1]-b[1]*c[0]);
}
int sign(std::int64_t value) { return (value > 0) - (value < 0); }
int reference(const std::array<IntegerPoint,5> &points)
{
  const auto &origin=points[0], &direction=points[1], &a=points[2], &b=points[3], &c=points[4];
  const auto ab=subtract(b,a), ac=subtract(c,a), ao=subtract(a,origin);
  const int denominator=sign(determinant(direction,ab,ac));
  const int numerator=sign(determinant(ab,ac,ao));
  if(denominator==0) { return numerator==0 ? 2 : 0; }
  for(int edge=0; edge<3; ++edge)
  {
    if(sign(determinant(direction,subtract(points[2+edge],origin),
      subtract(points[2+(edge+1)%3],origin)))==-denominator) { return 0; }
  }
  return numerator==0 || numerator==denominator ? 1 : 0;
}
}

TEST_CASE("predicates::exact_ray_triangle: integer_reference", "[igl][predicates]")
{
  // Verify classifications against independent bounded integer arithmetic.
  std::mt19937 generator(2549);
  std::uniform_int_distribution<int> sample(-4,4);
  for(int iteration=0; iteration<10000; ++iteration)
  {
    std::array<IntegerPoint,5> points;
    std::array<Eigen::RowVector3d,5> vectors;
    const double scale=std::ldexp(1.0, (iteration%7-3)*30);
    for(int i=0; i<5; ++i)
    {
      for(int axis=0; axis<3; ++axis)
      {
        points[i][axis]=sample(generator);
        vectors[i][axis]=double(points[i][axis])*scale;
      }
    }
    REQUIRE(int(igl::predicates::exact_ray_triangle_intersect(
      vectors[0],vectors[1],vectors[2],vectors[3],vectors[4]))==reference(points));
  }
}

TEST_CASE("predicates::exact_ray_triangle: integer_ordering", "[igl][predicates]")
{
  // Compare exact hit heights against an independent integer expression.
  std::mt19937 generator(916);
  std::uniform_int_distribution<int> sample(1,16);
  using Point=Eigen::RowVector3d;
  for(int iteration=0; iteration<5000; ++iteration)
  {
    const double scale=std::ldexp(1.0,(iteration%7-3)*30);
    const Point origin(0,0,0), direction(0,0,scale);
    const int radius1=sample(generator), radius2=sample(generator);
    const int z1=sample(generator), z2=sample(generator), z3=sample(generator);
    const int z4=sample(generator), z5=sample(generator), z6=sample(generator);
    const Point a1=Point(-radius1,-radius1,z1)*scale;
    Point b1=Point(radius1,-radius1,z2)*scale, c1=Point(0,radius1,z3)*scale;
    const Point a2=Point(-radius2,-radius2,z4)*scale;
    Point b2=Point(radius2,-radius2,z5)*scale, c2=Point(0,radius2,z6)*scale;
    if(iteration%2) { std::swap(b1,c1); }
    if(iteration%3) { std::swap(b2,c2); }
    const int expected=sign((z1+z2+2*z3)-(z4+z5+2*z6));
    REQUIRE(igl::predicates::exact_ray_triangle_hit_compare(
      origin,direction,a1,b1,c1,a2,b2,c2)==expected);
  }
}

TEST_CASE("predicates::exact_ray_triangle: exact_fallback", "[igl][predicates]")
{
  // Preserve a nonzero determinant lost by ordinary multiplication.
  using Point=Eigen::RowVector3d;
  using Result=igl::predicates::RayTriangleIntersection;
  const double epsilon=std::ldexp(1.0,-52);
  const Point origin(0,0,0), direction(0,0,1);
  const Point a(0,0,1), b(1,1+epsilon,1), c(1-epsilon,1,1);
  REQUIRE(igl::predicates::exact_ray_triangle_intersect(origin,direction,a,b,c)==Result::Hit);
  const Point a2(0,0,2), b2(1,1+epsilon,2), c2(1-epsilon,1,2);
  REQUIRE(igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a,b,c,a2,b2,c2)==-1);
  REQUIRE(igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a2,c2,b2,a,b,c)==1);
  REQUIRE(igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a,b,c,c,b,a)==0);

  // Keep a free direction even when adding it to the origin rounds it away.
  const double height=std::ldexp(1.0,60);
  const Point high_origin(0,0,height), high_a(-1,-1,height+256),
    high_b(1,-1,height+256), high_c(0,1,height+256);
  REQUIRE((high_origin+direction).eval()==high_origin);
  REQUIRE(igl::predicates::exact_ray_triangle_intersect(
    high_origin,direction,high_a,high_b,high_c)==Result::Hit);
}

TEST_CASE("predicates::exact_ray_triangle: float", "[igl][predicates]")
{
  // Promote float inputs exactly while using the double predicates backend.
  using Point=Eigen::RowVector3f;
  const Point origin(0,0,0), direction(0,0,1), a(-1,-1,1), b(1,-1,1), c(0,1,1);
  REQUIRE(igl::predicates::exact_ray_triangle_intersect(origin,direction,a,b,c)
    ==igl::predicates::RayTriangleIntersection::Hit);
  REQUIRE(igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a,b,c,c,b,a)==0);
}

TEST_CASE("predicates::exact_ray_triangle: filtered_boundary_handoff", "[igl][predicates]")
{
  // Preserve zero signs when any stage of the filtered classifier needs fallback.
  using Point=Eigen::RowVector3d;
  using Result=igl::predicates::RayTriangleIntersection;
  const Point direction(0,0,1);
  const std::array<Point,3> vertices{{Point(0,0,1),Point(1,0,1),Point(0,1,1)}};
  for(int rotation=0; rotation<3; ++rotation)
  {
    const Point a=vertices[rotation], b=vertices[(rotation+1)%3], c=vertices[(rotation+2)%3];
    for(const Point &surface : {a,Point((a+b)*0.5),Point((a+b+c)/3.0)})
    {
      const Point below=surface-direction, above=surface+direction;
      REQUIRE(igl::predicates::exact_ray_triangle_intersect(below,direction,a,b,c)==Result::Hit);
      REQUIRE(igl::predicates::exact_ray_triangle_intersect(surface,direction,a,b,c)==Result::Hit);
      REQUIRE(igl::predicates::exact_ray_triangle_intersect(above,direction,a,b,c)==Result::None);
      REQUIRE(igl::predicates::exact_ray_triangle_intersect(above,direction,c,b,a)==Result::None);
    }
  }
}

TEST_CASE("predicates::expansion: exact_capacity", "[igl][predicates]")
{
  // Exercise the public arithmetic using buffers with exactly the promised capacity.
  igl::predicates::exactinit();
  double large, tail;
  predicates_two_sum(1.0,std::ldexp(1.0,-54),&large,&tail);
  REQUIRE(large==1.0);
  REQUIRE(tail==std::ldexp(1.0,-54));
  predicates_two_diff(1.0,std::ldexp(1.0,-54),&large,&tail);
  REQUIRE(large==1.0);
  REQUIRE(tail==-std::ldexp(1.0,-54));
  const double epsilon=std::ldexp(1.0,-52);
  predicates_two_product(1.0+epsilon,1.0-epsilon,&large,&tail);
  REQUIRE(large==1.0);
  REQUIRE(tail==-epsilon*epsilon);
  const double left[1]={1.0}, right[1]={-1.0};
  double result[2];
  REQUIRE(predicates_expansion_sum(1,left,1,right,result)==1);
  REQUIRE(result[0]==0.0);
  const double small[1]={epsilon};
  REQUIRE(predicates_expansion_scale(1,small,epsilon,result)==1);
  REQUIRE(result[0]==epsilon*epsilon);
}
