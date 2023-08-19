#include <test_common.h>
#include <igl/exact_geodesic.h>

TEST_CASE("exact_geodesic: square", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(4,2);
  V << 0,0,
       1,0,
       1,1,
       0,1;
  Eigen::MatrixXi F(2,3);
  F << 0,1,2,
       0,2,3;
  Eigen::VectorXi VS(1);
  VS<<0;
  Eigen::VectorXi VT(4);
  VT<<0,1,2,3;
  Eigen::VectorXi FS,FT;
  Eigen::VectorXd D;
  igl::exact_geodesic(V,F,VS,FS,VT,FT,D);
  Eigen::VectorXd Dgt(4);
  Dgt<<0,1,1.4142135624,1;
  test_common::assert_near(D,Dgt,1e-10);
}
