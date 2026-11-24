#include "api.h"
#ifdef EXACT_RAYS_BASELINE
#include <igl/exact_ray_mesh_intersect.h>
#include <igl/exact_ray_triangle_intersect.h>
#include <igl/exact_ray_triangle_hit_compare.h>
namespace implementation = igl;
namespace baseline
#else
#include <igl/predicates/exact_ray_mesh_intersect.h>
#include <igl/predicates/exact_ray_triangle_intersect.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>
namespace implementation = igl::predicates;
namespace shewchuk
#endif
{
int classify(const Ray &ray, const Triangle &triangle)
{
  return int(implementation::exact_ray_triangle_intersect(
    ray.origin, ray.direction, triangle.a, triangle.b, triangle.c));
}
int compare(const Ray &ray, const Triangle &first, const Triangle &second)
{
  return implementation::exact_ray_triangle_hit_compare(ray.origin, ray.direction,
    first.a, first.b, first.c, second.a, second.b, second.c);
}
bool first(const Ray &ray, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces, int &face)
{
  return implementation::exact_ray_mesh_intersect(ray.origin,ray.direction,vertices,faces,face);
}
bool all(const Ray &ray, const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces, std::vector<int> &hits)
{
  return implementation::exact_ray_mesh_intersect_all(ray.origin,ray.direction,vertices,faces,hits);
}
double approximate(const Ray &ray, const Triangle &triangle)
{
  Eigen::MatrixXd vertices(3,3); vertices << triangle.a, triangle.b, triangle.c;
  Eigen::MatrixXi faces(1,3); faces << 0,1,2;
  igl::Hit<double> hit;
  implementation::exact_ray_mesh_intersect(ray.origin,ray.direction,vertices,faces,hit);
  return hit.t;
}
}
