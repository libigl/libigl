#include <test_common.h>
#include <igl/default_num_threads.h>
#include <igl/eytzinger_aabb.h>
#include <igl/box_simplices.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <igl/point_simplex_squared_distance.h>
#include <random>
#include <atomic>
#include <thread>

namespace
{
using Point=Eigen::RowVector3d;
namespace p=igl::predicates;
void require_feature(const Point &bary,int dimension,int id)
{
  for(int i=0;i<3;++i)
  {
    const bool active=dimension==2 || i==id || (dimension==1 && i==(id+1)%3);
    if(active) REQUIRE(bary(i)>0);
    else { REQUIRE(bary(i)==0); REQUIRE_FALSE(std::signbit(bary(i))); }
  }
}
struct Mesh
{
  Eigen::MatrixXd vertices,lower,upper;
  Eigen::MatrixXi faces,neighbors,neighbor_edges;
  Eigen::VectorXi leaves;
  p::MeshFeatureExposure exposure;
  void tree()
  {
    Eigen::MatrixXd minimum,maximum;
    igl::box_simplices(vertices,faces,minimum,maximum);
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
  Point reconstruct(Eigen::Index face,const Point &bary) const
  {
    Point result=Point::Zero();
    for(int i=0;i<3;++i) result+=bary(i)*vertices.row(faces(face,i));
    return result;
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
TEST_CASE("predicates::exact_point_mesh: cached_features_and_ties","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(8,3); mesh.faces.resize(8,3);
  mesh.vertices.topRows(4)<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
  mesh.faces.topRows(4)<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
  // Duplicate coordinates have different vertex ids, reversed edge directions
  // and different local feature ids. Tiny offsets also exercise exact ordering
  // when rounded squared distances cannot distinguish the two copies.
  for(int f=0;f<4;++f) for(int i=0;i<3;++i)
    mesh.faces(f+4,i)=mesh.faces(f,2-i)+4;
  const std::array<Point,10> queries{{Point(-1,-2,-3),Point(.25,-1,-2),
    Point(.25,.25,-1),Point(0,0,0),Point(.5,0,0),Point(.25,.25,0),
    Point(2,2,2),Point(.25,.25,.25),Point(0,0,std::ldexp(1.0,50)),
    Point(.25,.25,-std::ldexp(1.0,-52))}};
  for(const double offset:{0.0,std::ldexp(1.0,-52),.5})
  {
    mesh.vertices.bottomRows(4)=mesh.vertices.topRows(4);
    mesh.vertices.bottomRows(4).col(2).array()+=offset;
    mesh.tree();
    for(const auto &q:queries)
    {
      int best=0;
      for(int f=1;f<mesh.faces.rows();++f) if(mesh.compare(q,f,best)<0) best=f;
      REQUIRE(mesh.closest(q)==best);
    }
  }
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
  Eigen::VectorXd distances; Eigen::VectorXi triangles;
  Eigen::MatrixXd bary;
  p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,mesh.lower,
    mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
  REQUIRE(distances.size()==queries.rows()); REQUIRE(triangles.size()==8);
  REQUIRE(bary.rows()==8); REQUIRE(bary.cols()==3);
  for(int i=0;i<8;++i)
  {
    INFO("query "<<i);
    REQUIRE(distances(i)==Approx(expected[i]));
    REQUIRE(triangles(i)==0);
    require_feature(bary.row(i),dimensions[i],feature_ids[i]);
    REQUIRE(bary.row(i).sum()==Approx(1));
    Point reconstructed=Point::Zero();
    for(int j=0;j<3;++j) reconstructed+=bary(i,j)*mesh.vertices.row(mesh.faces(triangles(i),j));
    REQUIRE((queries.row(i)-reconstructed).norm()==Approx(std::abs(distances(i))));
    Eigen::Index triangle; Point weights;
    const Eigen::Vector3d query=queries.row(i).transpose();
    REQUIRE(p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,triangle,weights)==distances(i));
    REQUIRE(triangle==triangles(i)); REQUIRE(weights==bary.row(i));
    REQUIRE(p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,triangle,weights)==distances(i));
    REQUIRE(triangle==triangles(i)); REQUIRE(weights==bary.row(i));
    if(expected[i]==0) REQUIRE_FALSE(std::signbit(distances(i)));
  }
#ifndef IGL_STATIC_LIBRARY
  // Header-only float inputs and Eigen expressions preserve the same exact
  // classifications; all returned numerical constructions are doubles.
  Eigen::VectorXd float_distances; Eigen::VectorXi float_triangles;
  Eigen::MatrixXd float_bary;
  p::approximate_signed_distance_exact_sign(queries.cast<float>(),mesh.vertices.cast<float>(),
    mesh.faces,mesh.lower.cast<float>(),mesh.upper.cast<float>(),mesh.leaves,
    mesh.exposure,float_distances,float_triangles,float_bary);
  REQUIRE(float_distances==distances); REQUIRE(float_triangles==triangles);
  REQUIRE(float_bary==bary);
#endif
  queries(0,0)=std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary),std::invalid_argument);
  // Resizing reused output arrays must also handle empty batches.
  queries.resize(0,3);
  p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,mesh.lower,
    mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
  REQUIRE(distances.size()==0); REQUIRE(triangles.size()==0);
  REQUIRE(bary.rows()==0); REQUIRE(bary.cols()==3);
  queries.resize(1,2);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary),std::invalid_argument);
  queries.resize(1,3); queries.setZero(); mesh.exposure.edges.resize(0,3);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary),std::invalid_argument);
  mesh.vertices.resize(0,3); mesh.faces.resize(0,3); mesh.tree();
  mesh.exposure.vertices.resize(0);
  const Point zero=Point::Zero();
  Eigen::Index triangle; Point weights;
  REQUIRE_THROWS_AS(p::point_mesh_closest_point(zero,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,triangle,weights),std::invalid_argument);
  REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary),std::invalid_argument);
}

TEST_CASE("predicates::point_mesh_closest_point: cancellation_and_degeneracy","[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(3,3); mesh.faces.resize(1,3); mesh.faces<<0,1,2;
  const double thin=std::ldexp(1.,-40);
  mesh.vertices<<0,0,0, 1,1,0, 1,1+thin,0; mesh.tree();
  Eigen::Index triangle; Point bary;
  const auto evaluate=[&](const Point &query)
  { return p::point_mesh_closest_point(query,mesh.vertices,mesh.faces,mesh.lower,mesh.upper,mesh.leaves,triangle,bary); };
  const Point thin_query(.5,.5+thin/4,std::ldexp(1.,-20));
  double distance=evaluate(thin_query);
  require_feature(bary,2,0);
  REQUIRE(distance==std::ldexp(1.,-20));
  REQUIRE(bary==Point(.5,.25,.25));
  REQUIRE(mesh.reconstruct(triangle,bary)==Point(thin_query.x(),thin_query.y(),0));
  // The rounded constructed point has a different distance from the exact
  // point. Never use its norm to replace the selected feature's fraction.
  const double large=std::ldexp(1.,40),offset=std::ldexp(1.,-12);
  mesh.vertices<<0,0,0, 2*large,2*large,0, 0,0,1; mesh.tree();
  const Point query(large,large+offset,0);
  distance=evaluate(query);
  require_feature(bary,1,0);
  REQUIRE(distance==Approx(offset/std::sqrt(2.)).epsilon(1e-14));
  REQUIRE(std::abs((query-mesh.reconstruct(triangle,bary)).norm()-distance)>distance*.1);
  // Distance stays representable when the final squared-distance ratio
  // underflows. The polynomial intermediates themselves remain in range.
  const double scale=std::ldexp(1.,180),height=std::ldexp(1.,-600);
  mesh.vertices<<0,0,0, scale,0,0, 0,scale,0; mesh.tree();
  distance=evaluate(Point(scale/4,scale/4,height));
  require_feature(bary,2,0);
  REQUIRE(distance*distance==0); REQUIRE(distance==height);
  // Collinear and coincident triangles are valid for the unsigned wrapper.
  mesh.vertices<<0,0,0, 1,0,0, 2,0,0; mesh.tree();
  distance=evaluate(Point(.5,1,0));
  REQUIRE((bary.array()>0).count()==2);
  REQUIRE(distance==1); REQUIRE(mesh.reconstruct(triangle,bary)==Point(.5,0,0));
  mesh.vertices.setZero(); mesh.tree(); distance=evaluate(Point(1,2,2));
  require_feature(bary,0,0); REQUIRE(distance==3);
  REQUIRE(bary==Point(1,0,0));
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
    Eigen::Index triangle; Point bary;
    const double distance=p::approximate_signed_distance_exact_sign(query,mesh.vertices,
      mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,triangle,bary);
    REQUIRE((distance>0)-(distance<0)==-side); REQUIRE(triangle==0);
    require_feature(bary,2,0);
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
      Eigen::Index triangle; Point bary;
      const double distance=p::approximate_signed_distance_exact_sign(query,mesh.vertices,
        mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,triangle,bary);
      REQUIRE((distance>0)-(distance<0)==analytic_location(query,heights));
      Eigen::Index uncached_triangle; Point uncached_bary;
      REQUIRE(p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,uncached_triangle,uncached_bary)==distance);
      REQUIRE(uncached_triangle==triangle); REQUIRE(uncached_bary==bary);
      REQUIRE(triangle==mesh.closest(query));
      double reference=std::numeric_limits<double>::infinity();
      for(int face=0;face<mesh.faces.rows();++face)
      {
        double squared; Point point;
        igl::point_simplex_squared_distance<3>(query,mesh.vertices,mesh.faces,face,squared,point);
        reference=std::min(reference,squared);
      }
      REQUIRE(distance==Approx(analytic_location(query,heights)*std::sqrt(reference)).epsilon(1e-12));
      REQUIRE((query-mesh.reconstruct(triangle,bary)).norm()==Approx(std::abs(distance)).epsilon(1e-12));
    }
  }
}

TEST_CASE("predicates::approximate_signed_distance_exact_sign: parallel_batches","[igl][predicates]")
{
  const Mesh mesh=star({{-1,-1,-1,-1}});
  const std::array<Point,8> samples{{Point(0,0,.25),Point(0,0,-.25),Point(.25,.25,-.5),
    Point(1,0,-1),Point(1,1,-1),Point(0,0,-5),Point(0,-1,0),Point(0,0,-1)}};
  std::array<Eigen::Index,8> expected_triangle;
  std::array<Point,8> expected_bary;
  std::array<double,8> expected_distance;
  for(int i=0;i<8;++i) expected_distance[i]=p::approximate_signed_distance_exact_sign(samples[i],
    mesh.vertices,mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,
    expected_triangle[i],expected_bary[i]);
  Eigen::MatrixXd queries(10000,3);
  Eigen::VectorXd expected_distances(queries.rows());
  for(int i=0;i<queries.rows();++i)
  {
    queries.row(i)=samples[i%8]; expected_distances(i)=expected_distance[i%8];
  }
#ifndef IGL_STATIC_LIBRARY
  const auto caller=std::this_thread::get_id();
  std::atomic<bool> worker_read{false};
  // Observe input evaluation in a worker, so correctness alone cannot hide
  // an accidentally serial batch implementation.
  const auto observed_queries=Eigen::MatrixXd::NullaryExpr(queries.rows(),3,
    [&](Eigen::Index row,Eigen::Index col)
    {
      if(std::this_thread::get_id()!=caller) worker_read.store(true,std::memory_order_relaxed);
      return queries(row,col);
    });
  const auto serial_queries=observed_queries.topRows(9999);
#else
  // Static-library instantiations accept concrete matrices. Keep the same
  // batch correctness and exception checks; expression instrumentation is
  // available in the header-only build.
  const auto &observed_queries=queries;
  const Eigen::MatrixXd serial_queries=queries.topRows(9999);
#endif
  for(bool cached:{false,true})
  {
    Eigen::VectorXd distances; Eigen::VectorXi triangles;
    Eigen::MatrixXd bary;
    const auto evaluate=[&]
    {
      if(cached) p::approximate_signed_distance_exact_sign(observed_queries,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
      else p::approximate_signed_distance_exact_sign(observed_queries,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,distances,triangles,bary);
    };
#ifndef IGL_STATIC_LIBRARY
    worker_read=false;
#endif
    evaluate();
#ifndef IGL_STATIC_LIBRARY
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
    REQUIRE_FALSE(worker_read.load());
#else
    REQUIRE(worker_read.load()==(igl::default_num_threads()>1));
#endif
    worker_read=false;
#endif
    REQUIRE(distances==expected_distances);
    for(int i=0;i<queries.rows();i+=997)
    {
      REQUIRE(triangles(i)==expected_triangle[i%8]);
      REQUIRE(bary.row(i)==expected_bary[i%8]);
    }
    if(cached) p::approximate_signed_distance_exact_sign(serial_queries,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
    else p::approximate_signed_distance_exact_sign(serial_queries,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,distances,triangles,bary);
#ifndef IGL_STATIC_LIBRARY
    REQUIRE_FALSE(worker_read.load());
#endif
    REQUIRE(distances==expected_distances.head(9999));
    queries(9999,0)=std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(evaluate(),std::invalid_argument);
    queries.row(9999)=samples[9999%8];
  }
}

TEST_CASE("predicates::exact_point_mesh: symmetric_solid_tie", "[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(6,3);
  mesh.vertices.topRows(3)<<
    .09346180974327023,.31748807834095416,4.0811717497108745e-21,
    .79256025609371927,9.217454979427364e-21,-.43692309053096745,
    -1.9375255404802405e-21,-.74414186935724669,.5525544893106793;
  mesh.vertices.bottomRows(3)=-mesh.vertices.topRows(3);
  mesh.faces.resize(8,3);
  // The convex hull of three independent vectors and their negatives is a
  // closed octahedron. Orient each of its eight faces away from its center.
  for(int f=0;f<8;++f)
  {
    for(int c=0;c<3;++c) mesh.faces(f,c)=c+3*((f>>c)&1);
    const Point a=mesh.vertices.row(mesh.faces(f,0)),b=mesh.vertices.row(mesh.faces(f,1)),
      c=mesh.vertices.row(mesh.faces(f,2));
    if((b-a).cross(c-a).dot(a)<0) std::swap(mesh.faces(f,1),mesh.faces(f,2));
  }
  mesh.prepare();
  const Point query(0,0,0);
  // Opposite open faces are exactly equidistant. Tiny coordinates used to
  // lose low expansion bits and select face 7 instead of the lowest id 0.
  REQUIRE(mesh.compare(query,0,7)==0);
  REQUIRE(mesh.closest(query)==0);
  REQUIRE(p::exact_point_mesh_classify(query,mesh.vertices,mesh.faces,0,
    mesh.neighbors,mesh.neighbor_edges)==p::PointMeshLocation::INSIDE);
  REQUIRE(p::exact_point_mesh_classify(query,mesh.vertices,mesh.faces,0,
    mesh.exposure)==p::PointMeshLocation::INSIDE);
}

TEST_CASE("predicates::exact_point_mesh: reject_degenerate_signing_fold", "[igl][predicates]")
{
  for(bool degenerate_first:{false,true})
  {
    Mesh mesh;
    mesh.vertices.resize(5,3); mesh.vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,1, .5,0,0;
    mesh.faces.resize(6,3); mesh.faces<<0,1,4, 0,2,1, 0,4,3, 4,1,3, 0,3,2, 1,2,3;
    if(!degenerate_first) mesh.faces.row(0).swap(mesh.faces.row(5));
    mesh.tree();
    igl::triangle_triangle_adjacency(mesh.faces,mesh.neighbors,mesh.neighbor_edges);
    REQUIRE_THROWS_AS(p::exact_mesh_feature_exposure(mesh.vertices,mesh.faces,
      mesh.neighbors,mesh.neighbor_edges,mesh.exposure),std::invalid_argument);
    for(const Point &query:{Point(.25,-1,0),Point(.25,-1,-1),Point(.75,-.5,-.5),Point(.3,-2,-.1)})
    {
      Eigen::Index triangle; Point bary;
      REQUIRE_THROWS_AS(p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,triangle,bary),std::invalid_argument);
    }
  }
}

TEST_CASE("predicates::point_mesh_closest_point: rescaled_construction", "[igl][predicates]")
{
  for(int exponent:{-800,800})
  {
    const double scale=std::ldexp(1.,exponent);
    Mesh mesh;
    mesh.vertices=Eigen::MatrixXd::Zero(3,3);
    mesh.vertices(1,0)=scale; mesh.vertices(2,1)=scale;
    mesh.faces.resize(1,3); mesh.faces<<0,1,2;
    mesh.tree();
    const Point query=Point(.25,.25,1)*scale;
    Eigen::Index triangle; Point bary;
    const double distance=p::point_mesh_closest_point(query,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,triangle,bary);
    REQUIRE(triangle==0);
    require_feature(bary,2,0);
    REQUIRE(distance/scale==Approx(1.));
    const Point reconstructed=mesh.reconstruct(triangle,bary);
    REQUIRE(reconstructed(0)/scale==Approx(.25));
    REQUIRE(reconstructed(1)/scale==Approx(.25));
    REQUIRE(reconstructed(2)==0);
    REQUIRE(bary(0)==Approx(.5));
    if(exponent<0) REQUIRE(distance*distance==0);
    else REQUIRE(std::isinf(distance*distance));
  }
}

TEST_CASE("predicates::point_mesh_closest_point: barycentric_support_underflow", "[igl][predicates]")
{
  // The open-face weight at corner 1 is t^4/(M^4+M^2*t^2+t^4),
  // which rounds to zero. All polynomial intermediates stay in range.
  const double magnitude=std::ldexp(1.,140),tiny=std::ldexp(1.,-140);
  Mesh mesh;
  mesh.vertices.resize(3,3); mesh.vertices<<0,0,0, magnitude,tiny,0, 0,magnitude,tiny;
  mesh.faces.resize(1,3); mesh.faces<<0,1,2;
  mesh.tree();
  const Point query(0,tiny,0);
  Eigen::Index triangle; Point bary;
  const double distance=p::point_mesh_closest_point(query,mesh.vertices,mesh.faces,
    mesh.lower,mesh.upper,mesh.leaves,triangle,bary);
  REQUIRE(triangle==0);
  require_feature(bary,2,0);
  REQUIRE(bary(1)==std::numeric_limits<double>::denorm_min());
  REQUIRE(bary(2)/(tiny/magnitude)==Approx(1.));
  REQUIRE(distance/std::ldexp(1.,-420)==Approx(1.));
}

TEST_CASE("predicates::approximate_signed_distance_exact_sign: sign_underflow", "[igl][predicates]")
{
  const double scale=std::ldexp(1.,-1000),tilt=std::ldexp(1.,-40);
  for(int side:{-1,1})
  {
    Mesh mesh;
    mesh.vertices.resize(4,3);
    mesh.vertices<<0,0,0, scale,side*scale*tilt,0, 0,0,scale, 0,-scale,0;
    mesh.faces.resize(4,3); mesh.faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
    mesh.prepare();
    Eigen::MatrixXd queries(2,3); queries<<scale*tilt,0,scale/4, 0,0,0;
    for(bool cached:{false,true})
    {
      Eigen::VectorXd distances; Eigen::VectorXi triangles; Eigen::MatrixXd bary;
      if(cached) p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
      else p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,distances,triangles,bary);
      // The exact nonzero distance is about 2^-1080, below the rounding
      // midpoint to the smallest double. Its sign must remain observable.
      REQUIRE(distances(0)==-side*std::numeric_limits<double>::denorm_min());
      REQUIRE(triangles(0)==0); require_feature(bary.row(0),2,0);
      REQUIRE(distances(1)==0); REQUIRE_FALSE(std::signbit(distances(1)));
      require_feature(bary.row(1),0,0);
      Eigen::Index triangle; Point weights;
      const Point query=queries.row(0);
      const double distance=cached ? p::approximate_signed_distance_exact_sign(query,mesh.vertices,
        mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,triangle,weights) :
        p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
          mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,triangle,weights);
      REQUIRE(distance==distances(0)); REQUIRE(triangle==triangles(0));
      REQUIRE(weights==bary.row(0));
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
TEST_CASE("predicates::approximate_signed_distance_exact_sign: plain_object_outputs", "[igl][predicates]")
{
  Mesh mesh;
  mesh.vertices.resize(4,3); mesh.vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
  mesh.faces.resize(4,3); mesh.faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
  mesh.prepare();
  const double tiny=std::ldexp(1.,-160);
  Eigen::MatrixXd queries(4,3);
  queries<<.25,.25,tiny, .25,.25,-tiny, .25,.25,0, tiny,.25,-1;
  Eigen::Matrix<float,Eigen::Dynamic,1> distances;
  Eigen::Matrix<long long,Eigen::Dynamic,1> triangles;
  Eigen::Matrix<float,Eigen::Dynamic,3,Eigen::RowMajor> bary;
  for(bool cached:{false,true})
  {
    if(cached) p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,distances,triangles,bary);
    else p::approximate_signed_distance_exact_sign(queries,mesh.vertices,mesh.faces,
      mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,distances,triangles,bary);
    REQUIRE(distances.size()==4); REQUIRE(triangles.size()==4); REQUIRE(bary.rows()==4);
    REQUIRE(distances(0)==-std::numeric_limits<float>::denorm_min());
    REQUIRE(distances(1)==std::numeric_limits<float>::denorm_min());
    REQUIRE(distances(2)==0); REQUIRE_FALSE(std::signbit(distances(2)));
    REQUIRE(triangles(3)==0); REQUIRE((bary.row(3).array()>0).all());
    REQUIRE(bary(3,2)==std::numeric_limits<float>::denorm_min());
    const Point query=queries.row(3);
    Eigen::Index triangle;
    Eigen::Vector3f weights;
    const double distance=cached ? p::approximate_signed_distance_exact_sign(query,mesh.vertices,
      mesh.faces,mesh.lower,mesh.upper,mesh.leaves,mesh.exposure,triangle,weights) :
      p::approximate_signed_distance_exact_sign(query,mesh.vertices,mesh.faces,
        mesh.lower,mesh.upper,mesh.leaves,mesh.neighbors,mesh.neighbor_edges,triangle,weights);
    REQUIRE(distance==1); REQUIRE(triangle==0); REQUIRE(weights==bary.row(3).transpose());
    // A dynamic matrix can also hold the single-query row of barycentrics.
    Eigen::MatrixXd dynamic_bary;
    REQUIRE(p::point_mesh_closest_point(query,mesh.vertices,mesh.faces,mesh.lower,
      mesh.upper,mesh.leaves,triangle,dynamic_bary)==1);
    REQUIRE(dynamic_bary.rows()==1); REQUIRE(dynamic_bary.cols()==3);
    REQUIRE(dynamic_bary(0,2)==tiny);
  }
}
#endif
