#include <test_common.h>
#include <igl/eytzinger_aabb.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
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
