#include <test_common.h>
#include <igl/ray_triangle_intersection_parameters.h>

TEST_CASE("ray_triangle_intersection_parameters: basic", "[igl]")
{
  Eigen::RowVector3d A(0,0,0), B(1,0,0), C(0,1,0);
  Eigen::Vector3d source(0.25, 0.25, -2.0), dir(0,0,1);
  double t=0,u=0,v=0;
  REQUIRE(igl::ray_triangle_intersection_parameters(source, dir, A, B, C, t, u, v));
  REQUIRE(t == Approx(2.0));
  REQUIRE(u == Approx(0.25));
  REQUIRE(v == Approx(0.25));

  // Reconstruct the point two ways.
  const Eigen::RowVector3d p_ray = source.transpose() + t*dir.transpose();
  const Eigen::RowVector3d p_bary = A*(1-u-v) + B*u + C*v;
  REQUIRE((p_ray - p_bary).norm() == Approx(0.0).margin(1e-12));

  // No culling: a parameter outside the triangle is still returned.
  Eigen::Vector3d outside(2.0, 2.0, -1.0);
  REQUIRE(igl::ray_triangle_intersection_parameters(outside, dir, A, B, C, t, u, v));
  REQUIRE(u == Approx(2.0));
  REQUIRE(v == Approx(2.0));

  // Ray direction parallel to the triangle's (z=0) plane: no plane intersection.
  Eigen::Vector3d parallel(1,0,0);
  REQUIRE_FALSE(igl::ray_triangle_intersection_parameters(source, parallel, A, B, C, t, u, v));
}

TEST_CASE("ray_triangle_intersection_parameters: float", "[igl]")
{
  Eigen::Matrix<float,1,3> A(0,0,0), B(2,0,0), C(0,2,0);
  Eigen::Matrix<float,3,1> source(0.5f, 0.5f, 3.0f), dir(0,0,-1);
  float t=0,u=0,v=0;
  REQUIRE(igl::ray_triangle_intersection_parameters(source, dir, A, B, C, t, u, v));
  REQUIRE(t == Approx(3.0f));
  REQUIRE(u == Approx(0.25f));
  REQUIRE(v == Approx(0.25f));
}
