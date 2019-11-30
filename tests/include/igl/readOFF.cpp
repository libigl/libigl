#include <test_common.h>

TEST_CASE("readOFF: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOFF
    // directly...
    igl::read_triangle_mesh(test_common::data_path("cube.off"), V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 12);
    REQUIRE (F.cols() == 3);
}
