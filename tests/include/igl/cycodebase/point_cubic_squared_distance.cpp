#include <test_common.h>
#include <igl/cycodebase/point_cubic_squared_distance.h>

TEST_CASE("point_cubic_squared_distance: simple2d", "[igl/cycodebase]" )
{
  Eigen::MatrixXd C(4,2);
  C<<0,0,
    1,2,
    2,-2,
    3,0;
  Eigen::MatrixXd Q(3,2);
  //[1.5 0;2 0.5;2.5 1]
  Q<<
    1.5,0,
    2,0.5,
    2.5,1;
  Eigen::VectorXd sqrD,S;
  Eigen::MatrixXd K;
  igl::cycodebase::point_cubic_squared_distance(Q,C,sqrD,S,K);
  REQUIRE(sqrD.size() == 3);
  REQUIRE(S.size() == 3);
  REQUIRE(K.rows() == 3);
  REQUIRE(K.cols() == 2);
  // sqrd = [0;0.5;1.25]
  REQUIRE(sqrD[0] == Approx(0.0).margin(1e-12));
  REQUIRE(sqrD[1] == Approx(0.5).margin(1e-12));
  REQUIRE(sqrD[2] == Approx(1.25).margin(1e-12));
  // S = [0.5;0.5;1.0]
  REQUIRE(S[0] == Approx(0.5).margin(1e-12));
  REQUIRE(S[1] == Approx(0.5).margin(1e-12));
  REQUIRE(S[2] == Approx(1.0).margin(1e-12));
  // [1.5 0;1.5 0;3 0]
  REQUIRE(K(0,0) == Approx(1.5).margin(1e-12));
  REQUIRE(K(0,1) == Approx(0.0).margin(1e-12));
  REQUIRE(K(1,0) == Approx(1.5).margin(1e-12));
  REQUIRE(K(1,1) == Approx(0.0).margin(1e-12));
  REQUIRE(K(2,0) == Approx(3.0).margin(1e-12));
  REQUIRE(K(2,1) == Approx(0.0).margin(1e-12));
}
