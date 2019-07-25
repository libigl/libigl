#include <test_common.h>
#include <igl/path_to_edges.h>

TEST_CASE("igl_path_to_edges: basic_test", "[igl]")
{
  const Eigen::VectorXi I = (Eigen::VectorXi(6)<<0,1,2,3,4,5).finished();
  const Eigen::MatrixXi Eexpected = (Eigen::MatrixXi(4,5)<<0,1, 1,2, 2,3, 3,4, 4,5).finished();

  Eigen::MatrixXi Eactual;
  igl::path_to_edges(I, Eactual);

  test_common::assert_eq(Eactual, Eexpected);
}
