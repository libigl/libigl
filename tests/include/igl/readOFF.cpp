#include <test_common.h>
#include <igl/readOFF.h>

TEST_CASE("readOFF: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOFF(test_common::data_path("cube.off"), V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 12);
    REQUIRE (F.cols() == 3);
}
