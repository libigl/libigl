#include <test_common.h>
#include <igl/turning_number.h>
#include <igl/matlab_format.h>
#include <igl/PI.h>
#include <iostream>

TEST_CASE("turning_number: pentagon", "[igl]")
{
  Eigen::MatrixXd V(5, 2);
  V << 0, 0    ,
       1, 0    ,
       1, 1    ,
       0.5, 1.5,
       0, 1    ;
  const double t_ccw = igl::turning_number(V);
  // Reverse rows of V
  V = V.colwise().reverse().eval();
  const double t_cw = igl::turning_number(V);
  REQUIRE( abs(t_ccw - 1) < 1e-16 );
  REQUIRE( abs(t_cw - -1) < 1e-16 );
}

TEST_CASE("turning_number: heptagram", "[igl]")
{
  Eigen::MatrixXd V(7,2);
  for (int i = 0; i < 7; i++)
  {
      // Multiply by 3 to get the {7/3} skipping pattern
      double theta = 2.0 * igl::PI * (3 * i) / 7;
      V(i, 0) = std::cos(theta);
      V(i, 1) = std::sin(theta);
  }
  const double t_ccw = igl::turning_number(V);
  // Reverse rows of V
  V = V.colwise().reverse().eval();
  const double t_cw = igl::turning_number(V);
  REQUIRE( abs(t_ccw - 3) < 1e-16 );
  REQUIRE( abs(t_cw - -3) < 1e-16 );
}
