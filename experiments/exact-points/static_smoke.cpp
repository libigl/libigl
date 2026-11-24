#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
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
