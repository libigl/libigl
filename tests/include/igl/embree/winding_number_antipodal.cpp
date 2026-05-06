#include <test_common.h>
#include <igl/winding_number_antipodal.h>
#include <igl/WindingNumberAntipodalScene.h>
#include <igl/embree/EmbreeIntersector.h>
#include <igl/winding_number.h>
#include <igl/read_triangle_mesh.h>

#include <Eigen/Core>

#include <cmath>
#include <random>

TEST_CASE("winding_number_antipodal: closed mesh is integer-valued", "[igl/embree]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);

  // Build a grid of query points inside the AABB.
  const Eigen::RowVector3d mn = V.colwise().minCoeff();
  const Eigen::RowVector3d mx = V.colwise().maxCoeff();
  const Eigen::RowVector3d c  = 0.5 * (mn + mx);
  Eigen::MatrixXd P(8, 3);
  for (int i = 0; i < 8; ++i)
  {
    P(i,0) = c(0) + ((i & 1) ? 0.1 : -0.1) * (mx(0)-mn(0));
    P(i,1) = c(1) + ((i & 2) ? 0.1 : -0.1) * (mx(1)-mn(1));
    P(i,2) = c(2) + ((i & 4) ? 0.1 : -0.1) * (mx(2)-mn(2));
  }

  Eigen::VectorXd W_ref;
  igl::winding_number(V, F, P, W_ref);

  igl::embree::EmbreeIntersector e;
  e.init(V.cast<float>(), F.cast<int>());

  igl::WindingNumberAntipodalScene<double> scene(V, F);
  Eigen::VectorXd W;
  scene.winding_number(e, P, W);

  REQUIRE(W.size() == W_ref.size());
  for (int i = 0; i < W.size(); ++i)
  {
    REQUIRE(std::abs(W(i) - W_ref(i)) < 1e-3);
  }
}

TEST_CASE("winding_number_antipodal: single triangle has fractional contribution", "[igl/embree]")
{
  // Triangle in the z=0 plane, oriented CCW seen from +z.
  Eigen::MatrixXd V(3, 3);
  V << 0, 0, 0,
       1, 0, 0,
       0, 1, 0;
  Eigen::MatrixXi F(1, 3);
  F << 0, 1, 2;

  igl::embree::EmbreeIntersector e;
  e.init(V.cast<float>(), F.cast<int>());
  igl::WindingNumberAntipodalScene<double> scene(V, F);

  // Query points well above and below the triangle, near its centroid.
  const Eigen::RowVector3d c(1.0/3.0, 1.0/3.0, 0.0);
  Eigen::MatrixXd P(2, 3);
  P.row(0) = c + Eigen::RowVector3d(0, 0,  0.5);
  P.row(1) = c + Eigen::RowVector3d(0, 0, -0.5);

  Eigen::VectorXd W;
  scene.winding_number(e, P, W);

  // For a single oriented triangle, the GWN is roughly ±(solid angle / 4π).
  // For a point on the normal axis at height h above the centroid of an
  // equilateral-ish triangle, |W| stays well under 0.5 but is nonzero, and
  // the two sides have opposite signs.
  REQUIRE(std::abs(W(0)) > 1e-3);
  REQUIRE(std::abs(W(1)) > 1e-3);
  REQUIRE(W(0) * W(1) < 0.0);
}

TEST_CASE("winding_number_antipodal: random non-manifold soup matches reference", "[igl/embree]")
{
  // Fixed seed: deterministic across runs.
  std::mt19937 rng(0xA1B2C3D4u);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::uniform_int_distribution<int>     vidx(0, 9);

  // 10 random vertices.
  Eigen::MatrixXd V(10, 3);
  for (int i = 0; i < V.rows(); ++i)
    for (int d = 0; d < 3; ++d)
      V(i, d) = coord(rng);

  // 30 random triangles (degenerate ones, repeated indices, are skipped
  // and resampled so igl::winding_number's solid_angle stays well-defined).
  Eigen::MatrixXi F(30, 3);
  for (int t = 0; t < F.rows(); ++t)
  {
    int a, b, c;
    do { a = vidx(rng); b = vidx(rng); c = vidx(rng); }
    while (a == b || b == c || a == c);
    F(t, 0) = a; F(t, 1) = b; F(t, 2) = c;
  }

  // 100 random query points in a slightly larger box.
  Eigen::MatrixXd P(100, 3);
  std::uniform_real_distribution<double> qcoord(-1.5, 1.5);
  for (int i = 0; i < P.rows(); ++i)
    for (int d = 0; d < 3; ++d)
      P(i, d) = qcoord(rng);

  Eigen::VectorXd W_ref;
  igl::winding_number(V, F, P, W_ref);

  igl::embree::EmbreeIntersector e;
  e.init(V.cast<float>(), F.cast<int>());

  igl::WindingNumberAntipodalScene<double> scene(V, F);
  Eigen::VectorXd W;
  scene.winding_number(e, P, W);

  REQUIRE(W.size() == W_ref.size());
  for (int i = 0; i < W.size(); ++i)
  {
    REQUIRE(std::abs(W(i) - W_ref(i)) < 1e-5);
  }
}

TEST_CASE("winding_number_antipodal: free-function one-shot matches scene", "[igl/embree]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);

  const Eigen::RowVector3d c = 0.5 * (V.colwise().minCoeff() + V.colwise().maxCoeff());
  Eigen::MatrixXd P(4, 3);
  P.row(0) = c;
  P.row(1) = c + Eigen::RowVector3d(0.1, 0.0, 0.0);
  P.row(2) = c + Eigen::RowVector3d(0.0, 0.1, 0.0);
  P.row(3) = c + Eigen::RowVector3d(1e3, 0.0, 0.0);  // far outside

  igl::embree::EmbreeIntersector e;
  e.init(V.cast<float>(), F.cast<int>());

  Eigen::VectorXd W;
  igl::winding_number_antipodal(V, F, e, P, W);

  Eigen::VectorXd W_ref;
  igl::winding_number(V, F, P, W_ref);

  REQUIRE(W.size() == W_ref.size());
  for (int i = 0; i < W.size(); ++i)
    REQUIRE(std::abs(W(i) - W_ref(i)) < 1e-3);
}

TEST_CASE("winding_number_antipodal: float precision compiles", "[igl/embree]")
{
  using V3f = Eigen::Matrix<float, Eigen::Dynamic, 3>;
  using F3i = Eigen::Matrix<int,   Eigen::Dynamic, 3>;
  V3f V(3, 3);
  V << 0, 0, 0,
       1, 0, 0,
       0, 1, 0;
  F3i F(1, 3);
  F << 0, 1, 2;

  igl::embree::EmbreeIntersector e;
  e.init(V, F);
  igl::WindingNumberAntipodalScene<float> scene(V, F);

  Eigen::Matrix<float, Eigen::Dynamic, 3> P(1, 3);
  P << 1.0f/3.0f, 1.0f/3.0f, 0.5f;

  Eigen::VectorXf W;
  scene.winding_number(e, P, W);
  REQUIRE(std::abs(W(0)) > 1e-3f);
}
