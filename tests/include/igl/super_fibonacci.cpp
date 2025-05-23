#include <test_common.h>
#include <igl/super_fibonacci.h>

TEST_CASE("super_fibonacci: simple", "[igl]")
{
  Eigen::MatrixXd Q;
  igl::super_fibonacci(10,Q);
  REQUIRE(Q.rows() == 10);
  REQUIRE((Q.rowwise().norm()).eval().maxCoeff() == Approx(1.0));
  REQUIRE((Q.rowwise().norm()).eval().minCoeff() == Approx(1.0));
}
