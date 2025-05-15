#include <test_common.h>
#include <igl/winding_number.h>

TEST_CASE("winding_number: bunny", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("bunny.off"), V, F);
  Eigen::RowVectorXd V_mean = V.colwise().mean();
  Eigen::MatrixXd P = (V.rowwise() - V_mean) * 0.8;
  P.rowwise() += V_mean;
  Eigen::VectorXd W;
  igl::winding_number(V, F, P, W);
}
