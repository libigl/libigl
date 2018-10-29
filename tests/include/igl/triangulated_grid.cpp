#include <test_common.h>
#include <igl/triangulated_grid.h>
#include <igl/doublearea.h>

TEST_CASE("triangulated_grid: area", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  const int nx = 4;
  const int ny = 7;
  igl::triangulated_grid(nx,ny,V,F);
  REQUIRE (nx*ny == V.rows());
  REQUIRE (2*(nx-1)*(ny-1) == F.rows());
  Eigen::VectorXd dblA;
  igl::doublearea(V,F,dblA);
  REQUIRE (2.0 == Approx (dblA.array().sum()).margin(1e-10));
  const Eigen::VectorXd dblAgt = 
    Eigen::VectorXd::Constant(
      2*(nx-1)*(ny-1),
      1,
      2.0/double(2*(nx-1)*(ny-1)));
  test_common::assert_near(dblA,dblAgt,1e-10);
}
