#include <test_common.h>
#include <igl/is_edge_manifold.h>

TEST_CASE("is_edge_manifold: positive", "[igl]" "[slow]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path(param), V, F);
    REQUIRE ( igl::is_edge_manifold(F) );
  };

  test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}

TEST_CASE("is_edge_manifold: negative", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Known non-manifold mesh
  igl::read_triangle_mesh(test_common::data_path("truck.obj"), V, F);
  REQUIRE (! igl::is_edge_manifold(F) );
}
