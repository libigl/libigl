#include <test_common.h>
#include <igl/exact_geodesic.h>
#include <cmath>

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

TEST_CASE("exact_geodesic: disconnected components", "[igl]")
{
  Eigen::MatrixXd V(6,2);
  V << 0,0,
       1,0,
       0,1,
       2,0,
       3,0,
       2,1;
  Eigen::MatrixXi F(2,3);
  F << 0,1,2,
       3,4,5;
  Eigen::VectorXi VS(1);
  VS<<0;
  Eigen::VectorXi VT(2);
  VT<<1,5;
  Eigen::VectorXi FS,FT;
  Eigen::VectorXd D;
  igl::exact_geodesic(V,F,VS,FS,VT,FT,D);
  REQUIRE(D(0) == Approx(1.0));
  REQUIRE(std::isinf(D(1)));
}
