#include <test_common.h>
#include <igl/ray_mesh_intersect.h>

TEST_CASE("ray_mesh_intersect: one_triangle", "[igl]")
{
  Eigen::MatrixXd V(3,3);
  V.row(0) << 0.0, 0.0, 0.0;
  V.row(1) << 1.0, 0.0, 0.0;
  V.row(2) << 0.5, 1.0, 0.0;
  
  Eigen::MatrixXi F(1,3);
  F.row(0) << 0,1,2;
  
  Eigen::Vector3f source{0.5, 0.5, -1.0};
  Eigen::Vector3f direction{0.0, 0.0, 1.0};
  
  igl::Hit hit;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hit) == true);
  REQUIRE(hit.t == Approx(1.0));
  
  std::vector<igl::Hit> hits;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hits) == true);
  REQUIRE(hits.size() == 1);
  REQUIRE(hits.front().t == Approx(1.0));
}
