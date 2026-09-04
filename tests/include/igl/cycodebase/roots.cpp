#include <test_common.h>
#include <igl/cycodebase/roots.h>

TEST_CASE("roots: cubic", "[igl/cycodebase]" )
{
  // t^3 - 6*t^2 + 11*t - 6
  {
    Eigen::Vector<double,4> coef(4);
    coef << -6,11,-6,1;
    Eigen::Vector<double,3> R;
    const int nr = igl::cycodebase::roots(coef, 0.0, 4.0, R);
    REQUIRE(nr == 3);
    REQUIRE(R[0] == Approx(1.0).margin(1e-12));
    REQUIRE(R[1] == Approx(2.0).margin(1e-12));
    REQUIRE(R[2] == Approx(3.0).margin(1e-12));
  }

  // With dynamic size
  {
    Eigen::VectorXd coef(4);
    coef << -6,11,-6,1;
    {
      Eigen::VectorXd R;
      const int nr = igl::cycodebase::roots(coef, 0.0, 4.0, R);
      REQUIRE(nr == 3);
      REQUIRE(R.size() == 3);
      REQUIRE(R[0] == Approx(1.0).margin(1e-12));
      REQUIRE(R[1] == Approx(2.0).margin(1e-12));
      REQUIRE(R[2] == Approx(3.0).margin(1e-12));
    }
    // Only pick the first root
    {
      Eigen::VectorXd R;
      const int nr = igl::cycodebase::roots(coef, 0.0, 1.5, R);
      REQUIRE(nr == 1);
      REQUIRE(R.size() == 3);
      REQUIRE(R[0] == Approx(1.0).margin(1e-12));
      REQUIRE(std::isnan(R[1]));
      REQUIRE(std::isnan(R[2]));
    }
  }
}

