#include <test_common.h>
#include <igl/grid.h>
#include <igl/matlab_format.h>

TEST(grid,3d)
{
  Eigen::Vector3i res(3,3,3);
  Eigen::MatrixXd GV;
  igl::grid(res,GV);
  const Eigen::MatrixXd GVgt = 
    (Eigen::MatrixXd(27,3)<<
      0,  0,  0,
    0.5,  0,  0,
      1,  0,  0,
      0,0.5,  0,
    0.5,0.5,  0,
      1,0.5,  0,
      0,  1,  0,
    0.5,  1,  0,
      1,  1,  0,
      0,  0,0.5,
    0.5,  0,0.5,
      1,  0,0.5,
      0,0.5,0.5,
    0.5,0.5,0.5,
      1,0.5,0.5,
      0,  1,0.5,
    0.5,  1,0.5,
      1,  1,0.5,
      0,  0,  1,
    0.5,  0,  1,
      1,  0,  1,
      0,0.5,  1,
    0.5,0.5,  1,
      1,0.5,  1,
      0,  1,  1,
    0.5,  1,  1,
      1,  1,  1).finished();
  test_common::assert_eq(GV,GVgt);
}
