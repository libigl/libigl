#include <test_common.h>
#include <igl/cat.h>

#include <vector>

TEST_CASE("cat: rows", "[igl]")
{
  std::vector<Eigen::RowVector3i> rows = {
     Eigen::RowVector3i(1, 2, 3),
     Eigen::RowVector3i(4, 5, 6),
     Eigen::RowVector3i(7, 8, 9)
  };

  Eigen::MatrixXi actual;
  igl::cat(1,rows,actual);

  Eigen::Matrix3i expected;
  expected << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;  

  test_common::assert_eq(actual,expected);
}

TEST_CASE("cat: cols", "[igl]")
{
  std::vector<Eigen::Vector3i> cols = {
     Eigen::Vector3i(1, 4, 7),
     Eigen::Vector3i(2, 5, 8),
     Eigen::Vector3i(3, 6, 9)
  };

  Eigen::MatrixXi actual;
  igl::cat(2,cols,actual);

  Eigen::Matrix3i expected;
  expected << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;  

  test_common::assert_eq(actual,expected);
}

