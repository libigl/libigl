#include <test_common.h>

TEST(readOFF, simple) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOFF
    // directly...
    test_common::load_mesh("cube.off", V, F);
    ASSERT_EQ(8, V.rows());   //has 8 vertex
    ASSERT_EQ(3, V.cols());   //3D coordinates
    ASSERT_EQ(12, F.rows());  //has 6*2=12 facets
    ASSERT_EQ(3, F.cols());   //facets are triangles
}
