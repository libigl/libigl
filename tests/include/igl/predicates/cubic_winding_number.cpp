#include <test_common.h>
#include <igl/predicates/cubic_winding_number.h>
#include <iostream>

TEST_CASE("cubic_winding_number: simple", "[igl/predicates]")
{
  Eigen::Matrix<double,4,2> C;
  C<<0,0,
    1,1,
    2,-1,
    3,0;
  {
    Eigen::RowVector2d q(1.1,1.1);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(0.29147615882815 == Approx(w).epsilon(1e-12));
  }
  {
    Eigen::RowVector2d q(1.45,0.01);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(-0.502124394734 == Approx(w).epsilon(1e-12));
  }
  {
    Eigen::RowVector2d q(1.45,0);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(-0.5 == Approx(w).epsilon(1e-12));
  }
  // Stress test the numerics
  {
    Eigen::RowVector2d q(1.5+1e-15,0);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(0.5 == Approx(w).epsilon(1e-12));
  }
  {
    Eigen::RowVector2d q(1.5-1e-15,0);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(-0.5 == Approx(w).epsilon(1e-12));
  }
  {
    Eigen::RowVector2d q(1.5+1e-16,1e-16);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(0.5 == Approx(w).epsilon(1e-12));
  }
  {
    Eigen::RowVector2d q(1.5-1e-16,-1e-16);
    double w = igl::predicates::cubic_winding_number(C,q);
    REQUIRE(-0.5 == Approx(w).epsilon(1e-12));
  }
}

TEST_CASE("cubic_winding_number: degenerate", "[igl/predicates]")
{
  Eigen::RowVector2d q(1.2,1.0);
  {
    Eigen::Matrix<double,4,2> C;
    C<<0,0,
      0,0,
      0,0,
      0,0;
    REQUIRE( igl::predicates::cubic_winding_number(C,q) == 0.0);
  }
  {
    Eigen::Matrix<double,4,2> C;
    C<<0,0,
      1,0,
      2,0,
      3,0;
    REQUIRE(0.3087217355796 == Approx(igl::predicates::cubic_winding_number(C,q)).epsilon(1e-12));
  }
  {
    Eigen::Matrix<double,4,2> C;
    C<<0,0,
      1,1,
      3,3,
      9,9;
    REQUIRE(-0.4835565188700 == Approx(igl::predicates::cubic_winding_number(C,q)).epsilon(1e-12));
  }
}
