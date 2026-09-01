#include <test_common.h>
#include <igl/progressive_hulls.h>
#include <igl/block_self_intersections.h>
#include <igl/decimate_callback_types.h>
#include <igl/read_triangle_mesh.h>
#include <igl/centroid.h>
#include <vector>

TEST_CASE("progressive_hulls: reaches target face count on a closed mesh", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);
  const size_t target_m = F.rows()/2;

  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi J;
  const bool reached = igl::progressive_hulls(V,F,target_m,U,G,J);
  REQUIRE(reached);
  REQUIRE(static_cast<size_t>(G.rows()) == target_m);
  REQUIRE(U.allFinite());
}

TEST_CASE("progressive_hulls: output volume never shrinks below the input (outward placement)", "[igl]")
{
  // progressive hulls' entire point is that every collapse places the new
  // vertex outside the current mesh, so the enclosed volume should only grow
  // (up to floating-point/regularization noise), never shrink.
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);

  Eigen::Vector3d c0; double vol0;
  igl::centroid(V,F,c0,vol0);

  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi J;
  igl::progressive_hulls(V,F,F.rows()/2,U,G,J);

  Eigen::Vector3d c1; double vol1;
  igl::centroid(U,G,c1,vol1);

  REQUIRE(vol1 >= vol0 - 1e-6*std::abs(vol0));
}

TEST_CASE("progressive_hulls: matches deterministically across repeated runs", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);

  Eigen::MatrixXd U1,U2;
  Eigen::MatrixXi G1,G2;
  Eigen::VectorXi J1,J2;
  igl::progressive_hulls(V,F,F.rows()/3,U1,G1,J1);
  igl::progressive_hulls(V,F,F.rows()/3,U2,G2,J2);

  test_common::assert_near(U1,U2,0);
  REQUIRE(G1 == G2);
}

TEST_CASE("progressive_hulls: block_intersections flag matches decorator form", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"), V, F);
  const size_t target_m = F.rows()/2;

  // block_intersections=false is equivalent to the plain overload.
  {
    Eigen::MatrixXd U0,U1; Eigen::MatrixXi G0,G1; Eigen::VectorXi J0,J1;
    igl::progressive_hulls(V,F,target_m,U0,G0,J0);
    igl::progressive_hulls(V,F,target_m,/*block_intersections=*/false,U1,G1,J1);
    test_common::assert_near(U0,U1,0);
    REQUIRE(G0 == G1);
  }
  // block_intersections=true is equivalent to passing {block_self_intersections}.
  {
    std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> decs;
    decs.push_back(igl::block_self_intersections());
    Eigen::MatrixXd Ud,Ub; Eigen::MatrixXi Gd,Gb; Eigen::VectorXi Jd,Jb;
    igl::progressive_hulls(V,F,target_m,decs,Ud,Gd,Jd);
    igl::progressive_hulls(V,F,target_m,/*block_intersections=*/true,Ub,Gb,Jb);
    test_common::assert_near(Ud,Ub,0);
    REQUIRE(Gd == Gb);
  }
}
