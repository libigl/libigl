#include <igl/predicates/exact_ray_mesh_intersect.h>
#include <igl/predicates/exact_ray_triangle_intersect.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>
#include <stdexcept>

template<typename Point>
void check()
{
  const Point origin(0,0,0),direction(0,0,1);
  const Eigen::RowVector3d a(-1,-1,1),b(1,-1,1),c(0,1,1);
  Eigen::MatrixXd vertices(3,3); vertices<<a,b,c;
  Eigen::MatrixXi faces(1,3); faces<<0,1,2;
  int face=-1;
  std::vector<int> all_faces;
  igl::Hit<double> hit;
  std::vector<igl::Hit<double>> hits;
  const bool valid=igl::predicates::exact_ray_mesh_intersect(origin,direction,vertices,faces,face)
    && face==0 && igl::predicates::exact_ray_mesh_intersect_all(origin,direction,vertices,faces,all_faces)
    && all_faces==std::vector<int>{0}
    && igl::predicates::exact_ray_mesh_intersect(origin,direction,vertices,faces,hit) && hit.t==1
    && igl::predicates::exact_ray_mesh_intersect_all(origin,direction,vertices,faces,hits) && hits.size()==1
    && igl::predicates::exact_ray_triangle_intersect(origin,direction,a,b,c)==igl::predicates::RayTriangleIntersection::Hit
    && igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a,b,c,a,b,c)==0;
  if(!valid) { throw std::runtime_error("Static exact ray API smoke failed"); }
}
int main()
{
  check<Eigen::Vector3d>();
  check<Eigen::RowVector3d>();
  const Eigen::RowVector3f origin(0,0,0),direction(0,0,1),a(-1,-1,1),b(1,-1,1),c(0,1,1);
  return igl::predicates::exact_ray_triangle_intersect(origin,direction,a,b,c)!=igl::predicates::RayTriangleIntersection::Hit
    || igl::predicates::exact_ray_triangle_hit_compare(origin,direction,a,b,c,a,b,c)!=0;
}
