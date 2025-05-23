#include <test_common.h>
#include <igl/icosahedron.h>
#include <igl/doublearea.h>

TEST_CASE("super_fibonacci: simple", "[igl]")
{
  Eigen::MatrixXd Q;
  igl::super_fibonacci(10,Q);
  REQUIRE(Q.rows() == 10);
}
