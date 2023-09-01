#include <test_common.h>
#include <igl/edges_to_path.h>

TEST_CASE("edges_to_path: simple", "[igl]" )
{
  Eigen::MatrixXi E(3,2);
  E<<0,1,
     1,2,
     2,3;
  Eigen::VectorXi P,J,K;
  igl::edges_to_path(E,P,J,K);
  Eigen::VectorXi P_expected(4);
  P_expected<<0,1,2,3;
  test_common::assert_eq(P,P_expected);
}
