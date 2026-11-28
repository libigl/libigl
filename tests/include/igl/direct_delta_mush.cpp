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

TEST_CASE("direct_delta_mush: hinge with an unreferenced vertex", "[igl]")
{
  Eigen::MatrixXd V(5, 3), U, Omega;
  V << -2.064624, 10.402581, -0.046922,
        0.843939, 10.515798,  0.910529,
       -0.286211,  7.633783, -1.471316,
       -1.416073,  6.496789,  0.656846,
        1.,       -2.,        3.;
  Eigen::MatrixXi F(2, 3);
  F << 0, 3, 1,
       2, 0, 1;
  const Eigen::MatrixXd W = Eigen::MatrixXd::Ones(V.rows(), 1);
  const int p = GENERATE(1, 3);
  igl::direct_delta_mush_precomputation(V, F, W, p, 3., 1., 0.8, Omega);

  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> T(
    1, Eigen::Affine3d::Identity());
  const bool apply_transform = GENERATE(false, true);
  CAPTURE(p, apply_transform);
  CHECK((Omega.col(9) - W.col(0)).cwiseAbs().maxCoeff() < 1e-9);
  if (apply_transform)
  {
    T[0].translate(Eigen::Vector3d(1., -2., 3.));
    T[0].rotate(Eigen::AngleAxisd(0.6, Eigen::Vector3d::UnitZ()));
  }
  igl::direct_delta_mush(V, T, Omega, U);
  const Eigen::MatrixXd expected =
    (V * T[0].linear().transpose()).rowwise() + T[0].translation().transpose();
  test_common::assert_near(U, expected, 1e-9);
}
