#include <test_common.h>
#include <igl/cubic_split.h>

TEST_CASE("cubic_split: simple", "[igl]" )
{
  {
    Eigen::Matrix<double,4,2> C;
    C << 0,0,
         1,1,
         2,-1,
         3,0;
    double t = 0.5;
    Eigen::RowVector2d C01,C012,C0123,C123,C23;
    igl::cubic_split(C,t,C01,C012,C0123,C123,C23);
    // C01 = (0.5,0.5)
    // C012 = (1.0,0.25)
    // C0123 = (1.5,0.0)
    // C123 = (2.0,-0.25)
    // C23 = (2.5,-0.5)
    REQUIRE( C01.isApprox( Eigen::RowVector2d(0.5,0.5) ) );
    REQUIRE( C012.isApprox( Eigen::RowVector2d(1.0,0.25) ) );
    REQUIRE( C0123.isApprox( Eigen::RowVector2d(1.5,0.0) ) );
    REQUIRE( C123.isApprox( Eigen::RowVector2d(2.0,-0.25) ) );
    REQUIRE( C23.isApprox( Eigen::RowVector2d(2.5,-0.5) ) );
    Eigen::Matrix<double,4,2> C1,C2;
    igl::cubic_split(C,t,C1,C2);
    REQUIRE( C1.row(0).isApprox( C.row(0) ) );
    REQUIRE( C1.row(1).isApprox( C01 ) );
    REQUIRE( C1.row(2).isApprox( C012 ) );
    REQUIRE( C1.row(3).isApprox( C0123 ) );
    REQUIRE( C2.row(0).isApprox( C0123 ) );
    REQUIRE( C2.row(1).isApprox( C123 ) );
    REQUIRE( C2.row(2).isApprox( C23 ) );
    REQUIRE( C2.row(3).isApprox( C.row(3) ) );
  }
}
