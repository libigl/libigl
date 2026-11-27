#include <test_common.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_kernel.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/exactinit.h>
#include <array>
#include <cstdint>
#include <random>

namespace
{
using Integer=std::int64_t;
using IntegerPoint=std::array<Integer,3>;
using Point=Eigen::RowVector3d;
IntegerPoint subtract(const IntegerPoint &a,const IntegerPoint &b)
{ return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}}; }
Integer dot(const IntegerPoint &a,const IntegerPoint &b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
struct Reference { Integer numerator,denominator; int dimension,id; };
Reference reference(const IntegerPoint &q,const std::array<IntegerPoint,3> &p)
{
  Reference best{std::numeric_limits<Integer>::max(),1,0,0};
  const auto consider=[&](Integer n,Integer d,int dimension,int id)
  {
    if(best.numerator==std::numeric_limits<Integer>::max() ||
       n*best.denominator<best.numerator*d) best={n,d,dimension,id};
  };
  // Enumerate all feasible projections independently of Voronoi region tests.
  for(int i=0;i<3;++i)
  { const auto r=subtract(q,p[i]); consider(dot(r,r),1,0,i); }
  for(int i=0;i<3;++i)
  {
    const auto e=subtract(p[(i+1)%3],p[i]),r=subtract(q,p[i]);
    const Integer length=dot(e,e),t=dot(r,e);
    if(t>0 && t<length) consider(dot(r,r)*length-t*t,length,1,i);
  }
  const auto e=subtract(p[1],p[0]),f=subtract(p[2],p[0]),r=subtract(q,p[0]);
  const Integer x=dot(e,e),y=dot(f,f),z=dot(e,f),u=dot(r,e),v=dot(r,f);
  const Integer denominator=x*y-z*z,alpha=y*u-z*v,beta=x*v-z*u;
  if(denominator>0 && alpha>0 && beta>0 && alpha+beta<denominator)
    consider(dot(r,r)*denominator-(y*u*u-2*z*u*v+x*v*v),denominator,2,0);
  return best;
}
}
TEST_CASE("predicates::exact_point_triangle: integer_projection_oracle","[igl][predicates]")
{
  std::mt19937 generator(5917);
  std::uniform_int_distribution<int> coordinate(-4,4);
  for(int iteration=0;iteration<4000;++iteration)
  {
    IntegerPoint q;
    std::array<IntegerPoint,3> a,b;
    for(int axis=0;axis<3;++axis)
    {
      q[axis]=coordinate(generator);
      for(int i=0;i<3;++i) { a[i][axis]=coordinate(generator); b[i][axis]=coordinate(generator); }
    }
    const double scale=std::ldexp(1.0,(iteration%7-3)*20);
    const auto convert=[&](const IntegerPoint &p) { return Point(p[0]*scale,p[1]*scale,p[2]*scale); };
    const Point query=convert(q);
    const std::array<Point,3> aa{{convert(a[0]),convert(a[1]),convert(a[2])}},
      bb{{convert(b[0]),convert(b[1]),convert(b[2])}};
    const auto ar=reference(q,a),br=reference(q,b);
    int id;
    const auto feature=igl::predicates::exact_point_triangle_closest_feature(query,aa[0],aa[1],aa[2],id);
    REQUIRE(int(feature)==ar.dimension);
    REQUIRE(id==ar.id);
    namespace k=igl::predicates::exact_point_triangle_kernel;
    k::Fraction<k::Filter> distance;
    const int fused=k::filtered_feature_distance(k::point(query),
      {{k::point(aa[0]),k::point(aa[1]),k::point(aa[2])}},id,distance);
    if(fused>=0) { REQUIRE(fused==ar.dimension); REQUIRE(id==ar.id); }
    const Integer cross=ar.numerator*br.denominator-br.numerator*ar.denominator;
    const int expected=(cross>0)-(cross<0);
    REQUIRE(igl::predicates::exact_point_triangle_squared_distance_compare(query,
      aa[0],aa[1],aa[2],bb[0],bb[1],bb[2])==expected);
    REQUIRE(igl::predicates::exact_point_triangle_squared_distance_compare(query,
      bb[2],bb[1],bb[0],aa[2],aa[1],aa[0])==-expected);
  }
}
TEST_CASE("predicates::exact_point_triangle: boundaries_and_cancellation","[igl][predicates]")
{
  namespace p=igl::predicates;
  const Point a(0,0,0),b(1,0,0),c(0,1,0);
  const std::array<Point,7> queries{{Point(-1,-1,1),Point(2,0,1),Point(0,2,1),
    Point(.5,-1,1),Point(1,1,1),Point(-1,.5,1),Point(.25,.25,1)}};
  const int dimensions[7]={0,0,0,1,1,1,2},ids[7]={0,1,2,0,1,2,0};
  for(int i=0;i<7;++i)
  {
    int id;
    REQUIRE(int(p::exact_point_triangle_closest_feature(queries[i],a,b,c,id))==dimensions[i]);
    REQUIRE(id==ids[i]);
  }
  const double epsilon=std::ldexp(1.0,-52);
  const Point thin_b(1,1+epsilon,0),thin_c(1-epsilon,1,0),q(0,0,1);
  const Point upper_a(0,0,epsilon),upper_b=thin_b+upper_a,upper_c=thin_c+upper_a;
  REQUIRE(p::exact_point_triangle_squared_distance_compare(q,a,thin_b,thin_c,
    upper_a,upper_b,upper_c)==1);
  // The heights differ by less than one ulp of their rounded squared distances.
  const Point distant_query(0,0,std::ldexp(1.0,50)),near_a(0,0,epsilon),near_b=b+near_a,near_c=c+near_a;
  REQUIRE(p::exact_point_triangle_squared_distance_compare(distant_query,a,b,c,near_a,near_b,near_c)==1);
  REQUIRE(p::exact_point_triangle_squared_distance_compare(q,a,b,c,c,b,a)==0);
  int id;
  REQUIRE(p::exact_point_triangle_closest_feature(a,a,b,c,id)==p::ClosestPointFeature::VERTEX);
  REQUIRE(id==0);
  REQUIRE(p::exact_point_triangle_closest_feature(Point(.5,0,0),a,b,c,id)==p::ClosestPointFeature::EDGE);
  REQUIRE(id==0);
  const Eigen::Vector3f af(0,0,0),bf(1,0,0),cf(0,1,0),qf(.25f,.25f,1);
  REQUIRE(p::exact_point_triangle_closest_feature(qf,af,bf,cf,id)==p::ClosestPointFeature::FACE);
  REQUIRE(p::exact_point_triangle_squared_distance_compare(qf,af,bf,cf,cf,bf,af)==0);
}

TEST_CASE("predicates::exact_point_triangle: filter_certificates","[igl][predicates]")
{
  namespace k=igl::predicates::exact_point_triangle_kernel;
  igl::predicates::exactinit();
  std::mt19937 generator(81);
  std::uniform_real_distribution<double> mantissa(-1,1);
  std::uniform_int_distribution<int> exponent(-40,40);
  int certificates=0;
  for(int iteration=0;iteration<1000;++iteration)
  {
    k::Point q;
    std::array<k::Point,3> a,b;
    for(int axis=0;axis<3;++axis)
    {
      q[axis]=std::ldexp(mantissa(generator),exponent(generator));
      for(int i=0;i<3;++i)
      {
        a[i][axis]=std::ldexp(mantissa(generator),exponent(generator));
        b[i][axis]=std::ldexp(mantissa(generator),exponent(generator));
      }
    }
    if(iteration%2==0)
    {
      b=a;
      b[iteration%3][iteration%3]=std::nextafter(b[iteration%3][iteration%3],1.0);
    }
    int id,exact_id;
    const int feature=k::feature<k::Filter>(q,a,id);
    const int exact_feature=k::feature<k::Expansion>(q,a,exact_id);
    if(feature>=0) { REQUIRE(feature==exact_feature); REQUIRE(id==exact_id); }
    const int filtered=k::compare<k::Filter>(q,a,b),exact=k::compare<k::Expansion>(q,a,b);
    if(filtered!=2) { ++certificates; REQUIRE(filtered==exact); }
    REQUIRE(k::compare(q,a,b)==exact);
    int other_id;
    const int other_feature=k::feature<k::Expansion>(q,b,other_id);
    REQUIRE(k::compare(q,a,exact_feature,exact_id,b,other_feature,other_id)==exact);
  }
  REQUIRE(certificates>100);
}

TEST_CASE("predicates::exact_point_triangle: fused_distance_certificates","[igl][predicates]")
{
  namespace k=igl::predicates::exact_point_triangle_kernel;
  igl::predicates::exactinit();
  std::mt19937 generator(2151);
  std::uniform_real_distribution<double> mantissa(-1,1);
  std::uniform_int_distribution<int> exponent(-40,40),integer(-4,4);
  int certified=0,upper_certified=0,uncertain=0;
  for(int iteration=0;iteration<3000;++iteration)
  {
    k::Point q;
    std::array<k::Point,3> p;
    const auto coordinate=[&]()
    {
      return iteration<1500 ? std::ldexp(mantissa(generator),exponent(generator)) :
        std::ldexp(double(integer(generator)),(iteration%7-3)*20);
    };
    for(int axis=0;axis<3;++axis)
    {
      q[axis]=coordinate();
      for(int i=0;i<3;++i) p[i][axis]=coordinate();
    }
    if(iteration%13==0) q=p[0];
    if(iteration%17==0) p[2]=p[1];
    if(iteration%19==0) for(int axis=0;axis<3;++axis) p[2][axis]=(p[0][axis]+p[1][axis])/2;
    int id,exact_id;
    k::Fraction<k::Filter> fraction;
    int dimension=k::filtered_feature_distance(q,p,id,fraction);
    const int exact_dimension=k::feature<k::Expansion>(q,p,exact_id);
    if(dimension<0)
    {
      ++uncertain;
      dimension=k::feature(q,p,id);
      fraction=k::distance<k::Filter>(q,p,dimension,id);
    }
    else ++certified;
    REQUIRE(dimension==exact_dimension);
    REQUIRE(id==exact_id);
    const auto exact=k::distance<k::Expansion>(q,p,exact_dimension,exact_id);
    double upper;
    if(k::distance_upper_bound(fraction,upper))
    {
      ++upper_certified;
      REQUIRE(k::sign(k::Expansion(upper)*exact.denominator-exact.numerator)>=0);
    }
    auto other=p;
    double &perturbed=other[iteration%3][iteration%3];
    // nextafter(0,1) spans more than 1000 exponents relative to these meshes;
    // range rejection has its own test. Keep this certificate oracle in range.
    perturbed=perturbed==0 ? std::ldexp(1.,(iteration%7-3)*20-52) : std::nextafter(perturbed,1.0);
    int other_id;
    k::Fraction<k::Filter> other_distance;
    int other_dimension=k::filtered_feature_distance(q,other,other_id,other_distance);
    if(other_dimension<0)
    {
      other_dimension=k::feature(q,other,other_id);
      other_distance=k::distance<k::Filter>(q,other,other_dimension,other_id);
    }
    const int comparison=k::compare<k::Expansion>(q,p,dimension,id,other,other_dimension,other_id);
    REQUIRE(k::compare(q,p,dimension,id,fraction,other,other_dimension,other_id,other_distance)==comparison);
    // Certify N and D themselves, including the alternative edge formulas.
    for(const auto &pair:{std::make_pair(fraction.numerator,exact.numerator),
                          std::make_pair(fraction.denominator,exact.denominator)})
    {
      double error;
      if(!k::filter_error(pair.first,error)) continue;
      REQUIRE(k::sign(k::Expansion(k::up(pair.first.value+error))-pair.second)>=0);
      REQUIRE(k::sign(pair.second-k::Expansion(k::down(pair.first.value-error)))>=0);
    }
  }
  REQUIRE(certified>1000);
  REQUIRE(upper_certified>1000);
  REQUIRE(uncertain>100);
}

TEST_CASE("predicates::exact_point_triangle: distance_bound_refinement","[igl][predicates]")
{
  namespace k=igl::predicates::exact_point_triangle_kernel;
  igl::predicates::exactinit();
  const double epsilon=std::ldexp(1.0,-52);
  const k::Point q{{.5,.5+epsilon,.5}};
  const std::array<k::Point,3> p{{{{0,0,0}},{{2,2,0}},{{0,0,2}}}};
  int id;
  k::Fraction<k::Filter> fraction;
  REQUIRE(k::filtered_feature_distance(q,p,id,fraction)==2);
  double upper;
  REQUIRE_FALSE(k::distance_upper_bound(fraction,upper));
  // The existing interval/expansion path must resolve cancellation and zeros.
  REQUIRE(k::sign(k::distance<k::Expansion>(q,p,2,0).numerator)>0);
  REQUIRE(k::filtered_feature_distance(p[1],p,id,fraction)==-1);
  REQUIRE(k::feature(p[1],p,id)==0);
  REQUIRE(id==1);
  const auto zero=k::distance<k::Filter>(p[1],p,0,id);
  REQUIRE(k::distance_upper_bound(zero,upper));
  REQUIRE(upper==0);
  const double tiny=std::numeric_limits<double>::denorm_min();
  REQUIRE_FALSE(k::distance_upper_bound({k::Filter(tiny),k::Filter(1)},upper));
  REQUIRE_FALSE(k::distance_upper_bound({k::Filter(1),k::Filter(tiny)},upper));
  REQUIRE_FALSE(k::distance_upper_bound({k::Filter(std::numeric_limits<double>::infinity()),k::Filter(1)},upper));
}

TEST_CASE("predicates::exact_point_triangle: rescale_expansion_ties", "[igl][predicates]")
{
  namespace p=igl::predicates;
  const Point q(.0079102904640018235,-.39932183933392218,1.6866102000607417),
    a(-.48172845899923411,.3139446338798757,2.4651164674500633e-21),
    b(.87900262097467774,3.345215376473156e-21,-.18960187611408252),
    c(4.7334892242036771e-21,1.369257839175102,.56347102378932501);
  const Point aa(-.6317723228792238,.94272795638206652,-7.1722102170357587e-21),
    bb(.55475671348627009,8.0142140518424668e-21,-.95228725342804554),
    cc(4.8155278767464747e-21,-.90768105542844058,-.19858355535130057),zero(0,0,0);
  for(int exponent:{-400,0,400})
  {
    const double scale=std::ldexp(1.,exponent);
    const Point query=q*scale,x=a*scale,y=b*scale,z=c*scale;
    REQUIRE(p::exact_point_triangle_squared_distance_compare(query,x,y,z,z,y,x)==0);
    const Point first=aa*scale,second=bb*scale,third=cc*scale,
      opposite_first=-first,opposite_second=-second,opposite_third=-third;
    REQUIRE(p::exact_point_triangle_squared_distance_compare(zero,first,second,third,
      opposite_third,opposite_second,opposite_first)==0);
  }
  const Point x(1,0,0),y(0,1,0),z(std::ldexp(1.,-1070),0,1);
  REQUIRE_THROWS_AS(p::exact_point_triangle_squared_distance_compare(zero,x,y,z,z,y,x),std::range_error);
}

TEST_CASE("predicates::exact_point_triangle: filter_range_guard", "[igl][predicates]")
{
  namespace k=igl::predicates::exact_point_triangle_kernel;
  k::Point q{{.25,.25,1}};
  std::array<k::Point,3> triangle{{{{0,0,0}},{{1,0,0}},{{0,1,0}}}};
  REQUIRE(k::normal_filter_range(q,triangle));
  k::Fraction<k::Filter> fraction;
  int id;
  REQUIRE(k::filtered_feature_distance(q,triangle,id,fraction)==2);
  REQUIRE(fraction.numerator.value==1);
  for(auto &p:triangle) for(double &x:p) x=std::ldexp(x,-200);
  for(double &x:q) x=std::ldexp(x,-200);
  REQUIRE_FALSE(k::normal_filter_range(q,triangle));
  REQUIRE(k::filtered_feature_distance(q,triangle,id,fraction)==2);
  // The guarded path must flag the underflowed degree-six numerator. An
  // unchecked multiplication would silently return an exact-zero permanent.
  REQUIRE_FALSE(std::isfinite(fraction.numerator.permanent));
}
