#include <test_common.h>
#include <igl/cumprod.h>

TEST_CASE("cumprod: col_factorial", "[igl]")
{
  Eigen::Vector4d X(1,2,3,4);
  Eigen::Vector4d Y;
  igl::cumprod(X,1,Y);
  Eigen::Vector4d Ygt(1,2,6,24);
  test_common::assert_eq(Y,Ygt);
}

TEST_CASE("cumprod: row_factorial", "[igl]")
{
  Eigen::RowVector4d X(1,2,3,4);
  Eigen::RowVector4d Y;
  igl::cumprod(X,2,Y);
  Eigen::RowVector4d Ygt(1,2,6,24);
  test_common::assert_eq(Y,Ygt);
}
