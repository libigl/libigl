#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <stdexcept>

template<typename Point> void check()
{
  namespace p=igl::predicates;
  Eigen::MatrixXd vertices(4,3);
  vertices<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
  Eigen::MatrixXi faces(4,3),neighbors(4,3),edges(4,3);
  faces<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
  for(int f=0;f<4;++f) for(int i=0;i<3;++i)
    for(int g=0;g<4;++g) for(int j=0;j<3;++j)
      if(faces(f,i)==faces(g,(j+1)%3) && faces(f,(i+1)%3)==faces(g,j))
      { neighbors(f,i)=g; edges(f,i)=j; }
  p::MeshFeatureExposure exposure;
  p::exact_mesh_feature_exposure(vertices,faces,neighbors,edges,exposure);
  const Point query(.125,.125,.125);
  if(p::exact_point_mesh_classify(query,vertices,faces,0,neighbors,edges)!=p::PointMeshLocation::INSIDE ||
     p::exact_point_mesh_classify(query,vertices,faces,0,exposure)!=p::PointMeshLocation::INSIDE)
    throw std::runtime_error("Static local classification failed");
  Eigen::MatrixXd lower=Eigen::MatrixXd::Zero(7,3),upper=Eigen::MatrixXd::Ones(7,3);
  Eigen::VectorXi leaves(7); leaves<<-1,-1,-1,0,1,2,3;
  if(p::exact_point_mesh_closest_triangle(query,vertices,faces,lower,upper,leaves)!=0)
    throw std::runtime_error("Static closest triangle failed");
  const auto constructed=p::point_mesh_closest_point(query,vertices,faces,lower,upper,leaves);
  p::PointMeshClosestPoint signed_closest;
  p::PointMeshLocation location;
  const double distance=p::approximate_signed_distance_exact_sign(query,vertices,faces,
    lower,upper,leaves,exposure,signed_closest,location);
  Eigen::MatrixXd queries(1,3); queries.row(0)=Eigen::Map<const Eigen::RowVector3d>(query.data());
  Eigen::VectorXd distances; Eigen::VectorXi locations;
  std::vector<p::PointMeshClosestPoint> closest_points;
  p::approximate_signed_distance_exact_sign(queries,vertices,faces,lower,upper,leaves,
    exposure,distances,closest_points,locations);
  if(constructed.triangle!=0 || constructed.distance!=.125 || distance!=-.125 ||
     location!=p::PointMeshLocation::INSIDE || distances(0)!=distance ||
     locations(0)!=-1 || closest_points[0].triangle!=0)
    throw std::runtime_error("Static signed distance construction failed");
  if(p::approximate_signed_distance_exact_sign(query,vertices,faces,lower,upper,leaves,
    neighbors,edges,signed_closest,location)!=distance || location!=p::PointMeshLocation::INSIDE)
    throw std::runtime_error("Static uncached signed distance failed");
  p::approximate_signed_distance_exact_sign(queries,vertices,faces,lower,upper,leaves,
    neighbors,edges,distances,closest_points,locations);
  if(distances(0)!=distance || locations(0)!=-1 || closest_points[0].triangle!=0)
    throw std::runtime_error("Static uncached signed distance batch failed");
  const Eigen::RowVector3d a(0,0,0),b(1,0,0),c(0,1,0);
  int id;
  if(p::exact_point_triangle_closest_feature(query,a,b,c,id)!=p::ClosestPointFeature::FACE ||
     p::exact_point_triangle_squared_distance_compare(query,a,b,c,c,b,a)!=0)
    throw std::runtime_error("Static metric predicate failed");
}
int main()
{
  check<Eigen::RowVector3d>(); check<Eigen::Vector3d>();
  const Eigen::Vector3f query(.25,.25,1),a(0,0,0),b(1,0,0),c(0,1,0);
  int id;
  return igl::predicates::exact_point_triangle_closest_feature(query,a,b,c,id)!=igl::predicates::ClosestPointFeature::FACE ||
    igl::predicates::exact_point_triangle_squared_distance_compare(query,a,b,c,c,b,a)!=0;
}
