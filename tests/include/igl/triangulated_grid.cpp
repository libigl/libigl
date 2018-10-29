#include <test_common.h>
#include <igl/triangulated_grid.h>
#include <igl/doublearea.h>

TEST(triangulated_grid,area)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  const int nx = 4;
  const int ny = 7;
  igl::triangulated_grid(nx,ny,V,F);
  ASSERT_EQ(V.rows(),nx*ny);
  ASSERT_EQ(F.rows(),2*(nx-1)*(ny-1));
  Eigen::VectorXd dblA;
  igl::doublearea(V,F,dblA);
  ASSERT_NEAR(dblA.array().sum(),2.0,1e-10);
  const Eigen::VectorXd dblAgt = 
    Eigen::VectorXd::Constant(
      2*(nx-1)*(ny-1),
      1,
      2.0/double(2*(nx-1)*(ny-1)));
  test_common::assert_near(dblA,dblAgt,1e-10);
}
