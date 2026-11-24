#include <test_common.h>
#include <igl/cubic_is_flat.h>

TEST_CASE("cubic_is_flat: simple", "[igl]" )
{
  {
    Eigen::Matrix<double,4,2> C;
    C << 0,0,
         1,1,
         2,-1,
         3,0;
    REQUIRE( igl::cubic_is_flat(C,1e-2) == false);
  }
  {
    Eigen::Matrix<double,4,2> C;
    C << 0,0,
         1,0.1,
         2,-0.1,
         3,0;
    REQUIRE( igl::cubic_is_flat(C,1e-2) == true);
  }
  {
    Eigen::Matrix<double,4,2> C;
    C << 0,0,
         1,0.1,
         2,-0.1,
         3,0;
    REQUIRE( igl::cubic_is_flat(C,1e-4) == false);
  }
}

TEST_CASE("cubic_is_flat: degenerate", "[igl]" )
{
  Eigen::Matrix<double,4,2> C;
  C << 0,0,
    0,0,
    0,0,
    0,0;
  REQUIRE( igl::cubic_is_flat(C,1) == true);
}

