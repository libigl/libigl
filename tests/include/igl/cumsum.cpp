#include <test_common.h>
#include <igl/cumsum.h>

TEST_CASE("cumsum: col", "[igl]")
{
  Eigen::Vector4d X(1,2,3,4);
  Eigen::Vector4d Y;
  igl::cumsum(X,1,Y);
  Eigen::Vector4d Ygt(1,3,6,10);
  test_common::assert_eq(Y,Ygt);
}

TEST_CASE("cumsum: row", "[igl]")
{
  Eigen::RowVector4d X(1,2,3,4);
  Eigen::RowVector4d Y;
  igl::cumsum(X,2,Y);
  Eigen::RowVector4d Ygt(1,3,6,10);
  test_common::assert_eq(Y,Ygt);
}

