#include <test_common.h>
#include <igl/eytzinger_aabb.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <igl/point_simplex_squared_distance.h>
#include <random>

namespace
{
using Point=Eigen::RowVector3d;
namespace p=igl::predicates;
struct Mesh
{
  Eigen::MatrixXd vertices,lower,upper;
  Eigen::MatrixXi faces,neighbors,neighbor_edges;
  Eigen::VectorXi leaves;
  p::MeshFeatureExposure exposure;
  void tree()
  {
    Eigen::MatrixXd minimum(faces.rows(),3),maximum(faces.rows(),3);
    for(Eigen::Index f=0;f<faces.rows();++f)
    {
      minimum.row(f)=vertices.row(faces(f,0)); maximum.row(f)=minimum.row(f);
      for(int i=1;i<3;++i)
      {
        minimum.row(f)=minimum.row(f).cwiseMin(vertices.row(faces(f,i)));
        maximum.row(f)=maximum.row(f).cwiseMax(vertices.row(faces(f,i)));
      }
    }
    igl::eytzinger_aabb(minimum,maximum,lower,upper,leaves);
  }
  void prepare()
  {
    tree();
    igl::triangle_triangle_adjacency(faces,neighbors,neighbor_edges);
    p::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,exposure);
  }
  int compare(const Point &q,int a,int b) const
  {
    const Point a0=vertices.row(faces(a,0)),a1=vertices.row(faces(a,1)),a2=vertices.row(faces(a,2));
    const Point b0=vertices.row(faces(b,0)),b1=vertices.row(faces(b,1)),b2=vertices.row(faces(b,2));
    return p::exact_point_triangle_squared_distance_compare(q,a0,a1,a2,b0,b1,b2);
  }
  Eigen::Index closest(const Point &q) const
  { return p::exact_point_mesh_closest_triangle(q,vertices,faces,lower,upper,leaves); }
};
Mesh star(const std::array<double,4> &heights)
{
  Mesh mesh;
  mesh.vertices.resize(6,3);
  mesh.vertices<<0,0,0, 1,0,heights[0], 0,1,heights[1], -1,0,heights[2], 0,-1,heights[3], 0,0,-5;
  mesh.faces.resize(8,3);
  for(int i=0;i<4;++i)
  {
    mesh.faces.row(i)<<0,1+i,1+(i+1)%4;
    mesh.faces.row(i+4)<<5,1+(i+1)%4,1+i;
  }
  mesh.prepare(); return mesh;
}
int analytic_location(const Point &q,const std::array<double,4> &heights)
{
  const double radius=std::abs(q.x())+std::abs(q.y());
  if(radius>1) return 1;
  const double top=std::abs(q.x())*heights[q.x()>=0 ? 0 : 2]+
    std::abs(q.y())*heights[q.y()>=0 ? 1 : 3];
  const double bottom=top-5*(1-radius);
  if(q.z()>top || q.z()<bottom) return 1;
  if(q.z()==top || q.z()==bottom || radius==1) return 0;
  return -1;
}
}
TEST_CASE("predicates::exact_point_mesh: conservative_tree","[igl][predicates]")
{
  Mesh mesh;
  std::mt19937 generator(594);
  std::uniform_int_distribution<int> coordinate(-16,16);
  mesh.vertices.resize(75,3); mesh.faces.resize(25,3);
  for(int i=0;i<75;++i) for(int axis=0;axis<3;++axis) mesh.vertices(i,axis)=coordinate(generator)/8.0;
  for(int f=0;f<25;++f) mesh.faces.row(f)<<3*f,3*f+1,3*f+2;
  mesh.tree();
  for(int iteration=0;iteration<300;++iteration)
  {
    const Point q(coordinate(generator)/8.0,coordinate(generator)/8.0,coordinate(generator)/8.0);
    int best=0;
    for(int f=1;f<25;++f) if(mesh.compare(q,f,best)<0) best=f;
    REQUIRE(mesh.closest(q)==best);
  }
  // A rounded squared distance must not discard an exactly better candidate.
  mesh.vertices.resize(6,3); mesh.faces.resize(2,3);
  const double epsilon=std::ldexp(1.0,-52);
  mesh.vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,epsilon, 1,0,epsilon, 0,1,epsilon;
  mesh.faces<<0,1,2,3,4,5; mesh.tree();
  REQUIRE(mesh.closest(Point(0,0,std::ldexp(1.0,50)))==1);
  mesh.vertices.bottomRows(3)=mesh.vertices.topRows(3); mesh.tree();
  REQUIRE(mesh.closest(Point(.25,.25,1))==0);
  mesh.vertices.resize(0,3); mesh.faces.resize(0,3); mesh.tree();
  REQUIRE(mesh.closest(Point(0,0,0))==-1);
}
TEST_CASE("predicates::exact_point_mesh: local_and_cached_sign","[igl][predicates]")
{
  const std::array<std::array<double,4>,8> cases{{
    {{-1,-1,-1,-1}},{{1,1,1,1}},{{0,0,0,0}},{{1,-1,1,-1}},
    {{0,-1,0,-1}},{{0,1,0,1}},{{0,-1,0,0}},{{0,1,0,0}}}};
  const int labels[8]={1,-1,0,2,1,-1,1,-1};
  std::mt19937 generator(8991);
  std::uniform_int_distribution<int> coordinate(-12,12);
  for(int fixture=0;fixture<8;++fixture)
  {
    const Mesh mesh=star(cases[fixture]);
    REQUIRE(mesh.exposure.vertices(0)==labels[fixture]);
    // The supporting edge selected first by preprocessing is flat here;
    // the off-plane neighbor must supply the nonflat vertex's cached sign.
    if(fixture>=6) REQUIRE(mesh.exposure.edges(3,0)==0);
    for(int f=0;f<mesh.faces.rows();++f) for(int i=0;i<3;++i)
      REQUIRE(mesh.exposure.edges(f,i)==mesh.exposure.edges(mesh.neighbors(f,i),mesh.neighbor_edges(f,i)));
    for(int iteration=0;iteration<700;++iteration)
    {
      const Point q=iteration<6 ? Point(mesh.vertices.row(iteration)) :
        Point(coordinate(generator)/8.0,coordinate(generator)/8.0,coordinate(generator)/8.0);
      const auto closest=mesh.closest(q);
      const int expected=analytic_location(q,cases[fixture]);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,closest,mesh.neighbors,mesh.neighbor_edges))==expected);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,closest,mesh.exposure))==expected);
      if(iteration<20)
      {
        // All exact minimizers must give the same location, including ties.
        for(int face=0;face<8;++face) if(mesh.compare(q,face,closest)==0)
        {
          REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.neighbors,mesh.neighbor_edges))==expected);
          REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.exposure))==expected);
        }
      }
    }
    if(labels[fixture]!=2)
    {
      const Point q(0,0,labels[fixture]==-1 ? -.25 : .25);
      REQUIRE(mesh.closest(q)<4);
      int id;
      const int face=mesh.closest(q);
      const Point a=mesh.vertices.row(mesh.faces(face,0)),b=mesh.vertices.row(mesh.faces(face,1)),c=mesh.vertices.row(mesh.faces(face,2));
      REQUIRE(p::exact_point_triangle_closest_feature(q,a,b,c,id)==p::ClosestPointFeature::VERTEX);
      REQUIRE(mesh.faces(face,id)==0);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.exposure))==analytic_location(q,cases[fixture]));
    }
  }
}
TEST_CASE("predicates::exact_point_mesh: affine_and_orientation","[igl][predicates]")
{
  const std::array<double,4> heights{{0,-1,0,0}};
  std::mt19937 generator(130);
  std::uniform_int_distribution<int> coordinate(-12,12);
  for(int exponent:{-30,0,30})
  {
    Mesh mesh=star(heights);
    const double scale=std::ldexp(1.0,exponent);
    const auto transform=[&](const Point &q)
    { return Point(Point(-q.x()-2*q.y()+q.z()+4,q.y()+q.z()-2,q.z()+3)*scale); };
    for(int v=0;v<mesh.vertices.rows();++v) mesh.vertices.row(v)=transform(Point(mesh.vertices.row(v)));
    // The transform reverses orientation; reverse the faces to remain outward.
    mesh.faces.col(1).swap(mesh.faces.col(2)); mesh.prepare();
    REQUIRE(mesh.exposure.vertices(0)==1);
    for(int i=0;i<300;++i)
    {
      const Point original(coordinate(generator)/8.0,coordinate(generator)/8.0,coordinate(generator)/8.0);
      const Point q=transform(original);
      const auto face=mesh.closest(q);
      const int expected=analytic_location(original,heights);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.neighbors,mesh.neighbor_edges))==expected);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.exposure))==expected);
    }
  }
}
TEST_CASE("predicates::exact_point_mesh: irregular_stars_and_solid_angles","[igl][predicates]")
{
  std::mt19937 generator(517);
  std::uniform_real_distribution<double> sample(-1,1);
  for(int fixture=0;fixture<20;++fixture)
  {
    Mesh mesh;
    const int count=3+fixture%10;
    mesh.vertices.resize(count+2,3); mesh.faces.resize(2*count,3);
    mesh.vertices.row(0)<<0,0,0; mesh.vertices.row(count+1)<<0,0,-4;
    for(int i=0;i<count;++i)
    {
      const double angle=2*std::acos(-1.0)*i/count;
      mesh.vertices.row(i+1)<<std::cos(angle),std::sin(angle),sample(generator);
      mesh.faces.row(i)<<0,i+1,(i+1)%count+1;
      mesh.faces.row(i+count)<<count+1,(i+1)%count+1,i+1;
    }
    mesh.prepare();
    for(int iteration=0;iteration<100;++iteration)
    {
      const Point q(sample(generator)*1.5,sample(generator)*1.5,sample(generator)*2);
      long double total=0;
      for(int face=0;face<mesh.faces.rows();++face)
      {
        using LongPoint=Eigen::Matrix<long double,1,3>;
        const LongPoint a=mesh.vertices.row(mesh.faces(face,0)).cast<long double>()-q.cast<long double>();
        const LongPoint b=mesh.vertices.row(mesh.faces(face,1)).cast<long double>()-q.cast<long double>();
        const LongPoint c=mesh.vertices.row(mesh.faces(face,2)).cast<long double>()-q.cast<long double>();
        const long double determinant=a.x()*(b.y()*c.z()-b.z()*c.y())+
          a.y()*(b.z()*c.x()-b.x()*c.z())+a.z()*(b.x()*c.y()-b.y()*c.x());
        total+=2*std::atan2(determinant,a.norm()*b.norm()*c.norm()+a.dot(b)*c.norm()+b.dot(c)*a.norm()+c.dot(a)*b.norm());
      }
      const int expected=std::abs(total)>2*std::acos(-1.0L) ? -1 : 1;
      const auto face=mesh.closest(q);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.neighbors,mesh.neighbor_edges))==expected);
      REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.exposure))==expected);
    }
  }
}
TEST_CASE("predicates::exact_point_mesh: invalid_adjacency","[igl][predicates]")
{
  Mesh mesh=star({{-1,-1,-1,-1}});
  mesh.neighbors(0,0)=-1;
  REQUIRE_THROWS_AS(p::exact_mesh_feature_exposure(mesh.vertices,mesh.faces,mesh.neighbors,mesh.neighbor_edges,mesh.exposure),std::invalid_argument);
}

TEST_CASE("predicates::exact_point_mesh: hull_facet_boundary","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(8,3);
  mesh.vertices<<0,0,0, 2,0,0, 1,2,-1, -1,2,0,
    -2,0,-1, -1,-2,0, 1,-2,-1, 0,0,-5;
  mesh.faces.resize(12,3);
  for(int i=0;i<6;++i)
  {
    mesh.faces.row(i)<<0,i+1,(i+1)%6+1;
    mesh.faces.row(i+6)<<7,(i+1)%6+1,i+1;
  }
  mesh.prepare();
  // The center lies strictly inside a supporting hull facet, with a one-ray
  // normal cone. It is reachable despite not being a hull vertex or hull edge.
  REQUIRE(mesh.exposure.vertices(0)==int(p::FeatureExposure::EXPOSED_TO_OUTSIDE));
  const Point q(0,0,.25);
  const auto face=mesh.closest(q);
  REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.neighbors,mesh.neighbor_edges))==1);
  REQUIRE(int(p::exact_point_mesh_classify(q,mesh.vertices,mesh.faces,face,mesh.exposure))==1);
}

TEST_CASE("predicates::exact_point_mesh: cavity_shell","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(8,3); mesh.faces.resize(8,3);
  mesh.vertices<<0,0,0, 8,0,0, 0,8,0, 0,0,8,
    1,1,1, 2,1,1, 1,2,1, 1,1,2;
  mesh.faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3,
    4,5,6, 4,7,5, 4,6,7, 5,7,6;
  mesh.prepare();
  REQUIRE(mesh.exposure.vertices(0)==1);
  REQUIRE(mesh.exposure.vertices(4)==-1);
  const Point queries[5]={Point(1.125,1.125,1.125),Point(.5,.5,.5),Point(1,1,1),Point(1,1,.75),Point(0,0,-1)};
  const int expected[5]={1,-1,0,-1,1};
  for(int i=0;i<5;++i)
  {
    const auto face=mesh.closest(queries[i]);
    REQUIRE(int(p::exact_point_mesh_classify(queries[i],mesh.vertices,mesh.faces,face,mesh.neighbors,mesh.neighbor_edges))==expected[i]);
    REQUIRE(int(p::exact_point_mesh_classify(queries[i],mesh.vertices,mesh.faces,face,mesh.exposure))==expected[i]);
  }
}

TEST_CASE("predicates::approximate_signed_distance_exact_sign: batch_and_features","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(4,3); mesh.vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
  mesh.faces.resize(4,3); mesh.faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
  mesh.prepare();
  Eigen::MatrixXd queries(8,3);
  queries<<.25,.25,-1, .5,-1,-1, -1,-1,-1, .125,.125,.125,
    .25,.25,0, .5,0,0, 0,0,0, 2,0,0;
  const double expected[8]={1,std::sqrt(2.),std::sqrt(3.),-.125,0,0,0,1};
  const int dimensions[8]={2,1,0,2,2,1,0,0};
  const int feature_ids[8]={0,2,0,0,0,2,0,2};
  Eigen::VectorXd distances; Eigen::VectorXi locations;
  std::vector<p::PointMeshClosestPoint> closest;
  p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,mesh.lower,
    mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations);
  REQUIRE(distances.size()==queries.rows()); REQUIRE(closest.size()==8);
  for(int i=0;i<8;++i)
  {
    INFO("query "<<i);
    const auto &c=closest[i];
    REQUIRE(distances(i)==Approx(expected[i]));
    REQUIRE(locations(i)==(expected[i]>0)-(expected[i]<0));
    REQUIRE(c.triangle==0);
    REQUIRE(int(c.feature)==dimensions[i]); REQUIRE(c.feature_id==feature_ids[i]);
    REQUIRE(c.squared_distance==Approx(expected[i]*expected[i]));
    REQUIRE(c.distance==Approx(std::abs(expected[i])));
    REQUIRE(c.barycentric_coordinates.sum()==Approx(1));
    Point reconstructed=Point::Zero();
    for(int j=0;j<3;++j) reconstructed+=c.barycentric_coordinates(j)*mesh.vertices.row(mesh.faces(c.triangle,j));
    REQUIRE((c.closest_point-reconstructed).norm()<1e-14);
    REQUIRE((queries.row(i)-c.closest_point).norm()==Approx(c.distance));
    p::PointMeshClosestPoint single; p::PointMeshLocation location;
    const Eigen::Vector3d query=queries.row(i).transpose();
    REQUIRE(p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,single,location)==distances(i));
    REQUIRE(int(location)==locations(i)); REQUIRE(single.triangle==c.triangle);
    if(locations(i)==0) REQUIRE_FALSE(std::signbit(distances(i)));
  }
  // Header-only float inputs and Eigen expressions preserve the same exact
  // classifications; all returned numerical constructions are doubles.
  Eigen::VectorXd float_distances; Eigen::VectorXi float_locations;
  std::vector<p::PointMeshClosestPoint> float_closest;
  p::approximate_signed_distance_exact_sign(queries.cast<float>(),mesh.vertices.cast<float>(),
    mesh.faces,mesh.lower.cast<float>(),mesh.upper.cast<float>(),mesh.leaves,
    mesh.exposure,float_distances,float_closest,float_locations);
  REQUIRE(float_distances==distances); REQUIRE(float_locations==locations);
  queries(0,0)=std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations),std::invalid_argument);
  // Resizing reused output arrays must also handle empty batches.
  queries.resize(0,3);
  p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,mesh.lower,
    mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations);
  REQUIRE(distances.size()==0); REQUIRE(locations.size()==0); REQUIRE(closest.empty());
  queries.resize(1,2);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations),std::invalid_argument);
  queries.resize(1,3); queries.setZero(); mesh.exposure.edges.resize(0,3);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations),std::invalid_argument);
  mesh.vertices.resize(0,3); mesh.faces.resize(0,3); mesh.tree();
  mesh.exposure.vertices.resize(0);
  REQUIRE_THROWS_AS(p::point_mesh_closest_point(Point::Zero(),mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves),std::invalid_argument);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,closest,locations),std::invalid_argument);
}

TEST_CASE("predicates::point_mesh_closest_point: cancellation_and_degeneracy","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(3,3); mesh.faces.resize(1,3); mesh.faces<<0,1,2;
  const double thin=std::ldexp(1.,-40);
  mesh.vertices<<0,0,0, 1,1,0, 1,1+thin,0; mesh.tree();
  const auto evaluate=[&](const Point &query)
  { return p::point_mesh_closest_point(query,mesh.vertices,mesh.faces,mesh.lower,mesh.upper,mesh.leaves); };
  const Point thin_query(.5,.5+thin/4,std::ldexp(1.,-20));
  auto closest=evaluate(thin_query);
  REQUIRE(closest.feature==p::ClosestPointFeature::FACE);
  REQUIRE(closest.distance==std::ldexp(1.,-20));
  REQUIRE(closest.barycentric_coordinates==Point(.5,.25,.25));
  REQUIRE(closest.closest_point==Point(thin_query.x(),thin_query.y(),0));
  // The rounded constructed point has a different distance from the exact
  // point. Never use its norm to replace the selected feature's fraction.
  const double large=std::ldexp(1.,40),offset=std::ldexp(1.,-12);
  mesh.vertices<<0,0,0, 2*large,2*large,0, 0,0,1; mesh.tree();
  const Point query(large,large+offset,0);
  closest=evaluate(query);
  REQUIRE(closest.feature==p::ClosestPointFeature::EDGE);
  REQUIRE(closest.feature_id==0);
  REQUIRE(closest.distance==Approx(offset/std::sqrt(2.)).epsilon(1e-14));
  REQUIRE(closest.squared_distance==offset*offset/2);
  REQUIRE(std::abs((query-closest.closest_point).norm()-closest.distance)>closest.distance*.1);
  // Distance stays representable when the final squared-distance ratio
  // underflows. The polynomial intermediates themselves remain in range.
  const double scale=std::ldexp(1.,180),height=std::ldexp(1.,-600);
  mesh.vertices<<0,0,0, scale,0,0, 0,scale,0; mesh.tree();
  closest=evaluate(Point(scale/4,scale/4,height));
  REQUIRE(closest.feature==p::ClosestPointFeature::FACE);
  REQUIRE(closest.squared_distance==0); REQUIRE(closest.distance==height);
  // Collinear and coincident triangles are valid for the unsigned wrapper.
  mesh.vertices<<0,0,0, 1,0,0, 2,0,0; mesh.tree();
  closest=evaluate(Point(.5,1,0));
  REQUIRE(closest.feature==p::ClosestPointFeature::EDGE);
  REQUIRE(closest.distance==1); REQUIRE(closest.closest_point==Point(.5,0,0));
  mesh.vertices.setZero(); mesh.tree(); closest=evaluate(Point(1,2,2));
  REQUIRE(closest.feature==p::ClosestPointFeature::VERTEX);
  REQUIRE(closest.feature_id==0); REQUIRE(closest.distance==3);
  REQUIRE(closest.barycentric_coordinates==Point(1,0,0));
  REQUIRE_THROWS_AS(evaluate(Point(0,0,std::numeric_limits<double>::infinity())),std::invalid_argument);
}

TEST_CASE("predicates::approximate_signed_distance_exact_sign: near_surface","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(4,3); mesh.vertices<<0,0,0, 2,2,0, 0,0,2, 0,1,0;
  mesh.faces.resize(4,3); mesh.faces<<0,1,2, 0,3,1, 0,2,3, 1,3,2;
  mesh.prepare();
  const double offset=std::ldexp(1.,-52);
  for(int side:{-1,0,1})
  {
    const Point query(.5,.5+side*offset,.5);
    p::PointMeshClosestPoint closest; p::PointMeshLocation location;
    const double distance=p::approximate_signed_distance_exact_sign(query,mesh.vertices,
      mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,closest,location);
    REQUIRE(int(location)==-side); REQUIRE(closest.triangle==0);
    REQUIRE(closest.feature==p::ClosestPointFeature::FACE);
    if(side==0) { REQUIRE(distance==0); REQUIRE_FALSE(std::signbit(distance)); }
    else REQUIRE(distance/offset==Approx(-side/std::sqrt(2.)).epsilon(1e-14));
  }
}

TEST_CASE("predicates::approximate_signed_distance_exact_sign: analytic_stars","[igl][predicates]")
{
  std::mt19937 generator(2945);
  std::uniform_real_distribution<double> sample(-1.5,1.5);
  for(const std::array<double,4> heights : {std::array<double,4>{{-1,-1,-1,-1}},
    {{1,1,1,1}},{{0,0,0,0}},{{1,-1,1,-1}}})
  {
    const Mesh mesh=star(heights);
    for(int iteration=0;iteration<100;++iteration)
    {
      const Point query(sample(generator),sample(generator),sample(generator));
      p::PointMeshClosestPoint closest; p::PointMeshLocation location;
      const double distance=p::approximate_signed_distance_exact_sign(query,mesh.vertices,
        mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,closest,location);
      REQUIRE(int(location)==analytic_location(query,heights));
      REQUIRE(closest.triangle==mesh.closest(query));
      double reference=std::numeric_limits<double>::infinity();
      for(int face=0;face<mesh.faces.rows();++face)
      {
        double squared; Point point;
        igl::point_simplex_squared_distance<3>(query,mesh.vertices,mesh.faces,face,squared,point);
        reference=std::min(reference,squared);
      }
      REQUIRE(closest.squared_distance==Approx(reference).epsilon(1e-12));
      REQUIRE(distance==Approx(int(location)*std::sqrt(reference)).epsilon(1e-12));
      REQUIRE((query-closest.closest_point).norm()==Approx(closest.distance).epsilon(1e-12));
    }
  }
}
