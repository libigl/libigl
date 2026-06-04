#include <test_common.h>
#include <igl/cycodebase/box_cubic.h>

TEST_CASE("box_cubic: simple", "[igl/cycodebase]" )
{
  Eigen::MatrixXd C(4,2);
  C<<0,0,
    1,2,
    2,-2,
    3,0;
  Eigen::RowVectorXd B1,B2;
  igl::cycodebase::box_cubic(C,B1,B2);
  REQUIRE(B1.size() == 2);
  REQUIRE(B2.size() == 2);
  REQUIRE(B1[0] == Approx(0.0).margin(1e-12));
  REQUIRE(B1[1] == Approx(-0.57735026918962584).margin(1e-12));
  REQUIRE(B2[0] == Approx(3.0).margin(1e-12));
  REQUIRE(B2[1] == Approx(0.57735026918962584).margin(1e-12));
}

