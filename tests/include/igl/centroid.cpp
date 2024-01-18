#include <test_common.h>
#include <igl/centroid.h>
#include <iostream>

TEST_CASE("centroid: ", "[igl]" )
{
  Eigen::MatrixXd V(4,3);
  V<<0,0,0,
    1,0,0,
    0,1,0,
    0,0,1;
  Eigen::MatrixXi F(4,3);
    F<< 
    0,1,3,
    0,2,1,
    0,3,2,
    1,2,3;
  {
    Eigen::RowVector3d c;
    igl::centroid(V,F,c);
    REQUIRE(c(0) == Approx(0.25));
    REQUIRE(c(1) == Approx(0.25));
    REQUIRE(c(2) == Approx(0.25));
  }
  {
    Eigen::RowVectorXd c;
    igl::centroid(V,F,c);
    REQUIRE(c(0) == Approx(0.25));
    REQUIRE(c(1) == Approx(0.25));
    REQUIRE(c(2) == Approx(0.25));
  }
}

