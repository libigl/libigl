#include <test_common.h>

TEST(readOBJ, simple) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("cube.obj", V, F);
    ASSERT_EQ(8, V.rows());
    ASSERT_EQ(12, F.rows());
}
