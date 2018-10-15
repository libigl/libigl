#include <test_common.h>

TEST(readDMAT, Comp) {
    Eigen::MatrixXd N1, N2;
    test_common::load_matrix("duplicated_faces_N1.dmat", N1);
    test_common::load_matrix("duplicated_faces_N2.dmat", N2);

    ASSERT_EQ(N1.rows(), N2.rows());
    ASSERT_EQ(N1.cols(), N2.cols());
    ASSERT_FALSE(((N1-N2).array() != 0.0).all());

    const size_t rows = N1.rows();
    const size_t cols = N1.cols();
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            ASSERT_FLOAT_EQ(N1(i,j), N2(i,j));
        }
    }
}
