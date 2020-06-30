#include <test_common.h>
#include <igl/direct_delta_mush.h>
#include <igl/adjacency_list.h>
#include <iostream>

TEST_CASE("direct_delta_mush: cube", "[igl]")
{
  Eigen::MatrixXd V, W, U, Omega;
  Eigen::MatrixXi F;

  // Test that direct delta mush with identity transform reproduces the rest state of the mesh.
  // For simplicity, testing on a cube of dimensions 1.0 x 1.0 x 1.0,
  // with all vertices bound strictly to one and only one imaginary bone (weight is a column of 1s)
  igl::read_triangle_mesh(test_common::data_path("cube.off"), V, F);
  W = Eigen::MatrixXd::Ones(V.rows(), 1);

  // Parameters such as p, lambda, kappa, alpha do not matter given identity transform
  igl::direct_delta_mush_precomputation(V, F, W, 1, 1., 0.5, 0.5, Omega);

  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> T_list;
  T_list.push_back(Eigen::Affine3d::Identity());
  igl::direct_delta_mush(V, T_list, Omega, U);

  test_common::assert_near(U, V, 1e-4);
}
