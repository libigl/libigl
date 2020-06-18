#include <test_common.h>
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>
#include <iostream>

TEST_CASE("dijkstra: cube", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.off"), V, F);

  std::vector<std::vector<int>> VV;
  igl::adjacency_list(F, VV);

  Eigen::VectorXd min_distance;
  Eigen::VectorXi previous;
  igl::dijkstra(V, VV, 0, {7}, min_distance, previous);

  REQUIRE(min_distance(0) == 0);
  REQUIRE(min_distance(7) == Approx(sqrt(2)).margin(1e-10));
}

TEST_CASE("dijkstra: discrete distances", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.off"), V, F);

  std::vector<std::vector<int>> VV;
  igl::adjacency_list(F, VV);

  Eigen::VectorXi min_distance, previous;
  int out = igl::dijkstra(0, {3}, VV, min_distance, previous);
  REQUIRE(out != -1);
  REQUIRE(out == 3);
  REQUIRE(min_distance[3] == 1);
  REQUIRE(min_distance[0] == 0);
}

