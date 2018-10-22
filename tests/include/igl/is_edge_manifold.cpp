#include <test_common.h>
#include <igl/is_edge_manifold.h>


TEST_CASE("is_edge_manifold: positive", "[igl]")
{
	const auto test_case = [](const std::string &param)
	{
	  Eigen::MatrixXd V;
	  Eigen::MatrixXi F;
	  test_common::load_mesh(param, V, F);
	  REQUIRE ( igl::is_edge_manifold(F) );
	};

	test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}

TEST_CASE("is_edge_manifold: negative", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Known non-manifold mesh
  test_common::load_mesh("truck.obj", V, F);
  REQUIRE (! igl::is_edge_manifold(F) );
}
