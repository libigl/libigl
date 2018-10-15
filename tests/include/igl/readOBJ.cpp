#include <test_common.h>

TEST(readOBJ, simple) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOBJ
    // directly...
    test_common::load_mesh("cube.obj", V, F);
    ASSERT_EQ(8, V.rows());
    ASSERT_EQ(12, F.rows());
}
