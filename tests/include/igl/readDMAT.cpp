#include <test_common.h>

TEST_CASE("readDMAT: Comp", "[igl]")
{
    Eigen::MatrixXd N1, N2;
    test_common::load_matrix("duplicated_faces_N1.dmat", N1);
    test_common::load_matrix("duplicated_faces_N2.dmat", N2);

    REQUIRE (N2.rows() == N1.rows());
    REQUIRE (N2.cols() == N1.cols());
    REQUIRE (!((N1-N2).array() != 0.0).all());

    const size_t rows = N1.rows();
    const size_t cols = N1.cols();
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            REQUIRE (N2(i,j) == Approx(N1(i,j)));
        }
    }
}
