#include <test_common.h>

TEST_CASE("readOBJ: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOBJ
    // directly...
    test_common::load_mesh("cube.obj", V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (F.rows() == 12);
}
