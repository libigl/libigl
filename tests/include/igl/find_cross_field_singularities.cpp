#include <test_common.h>
#include <igl/find_cross_field_singularities.h>
#include <Eigen/Core>

TEST_CASE("find_cross_field_singularities", "[igl]")
{
    constexpr size_t N=8;

    Eigen::MatrixXd V(N + 2, 3);
    Eigen::MatrixXi F(N * 2, 3);
    Eigen::MatrixXd PD1(N * 2, 3);
    Eigen::MatrixXd PD2(N * 2, 3);

    V.row(0) << 0, 0, 1;
    V.row(1) << 0, 0, -1;
    for (size_t i=0; i<N; i++) {
        V.row(i + 2) << std::cos(2 * M_PI * i / N), std::sin(2 * M_PI * i / N), 0;
        F.row(i) << 0, i + 2, (i + 1) % N + 2;
        F.row(i + N) << 1, (i + 1) % N + 2, i + 2;
    }

    for (size_t i=0; i<N; i++) {
        PD1.row(i) = V.row(F(i, 2)) - V.row(F(i, 1));
        PD1.row(i + N) = PD1.row(i);
        PD2.row(i) = V.row(0) - (V.row(F(i, 1)) + V.row(F(i, 2))) / 2;
        PD2.row(i + N) = V.row(1) - (V.row(F(i + N, 1)) + V.row(F(i + N, 2))) / 2;
    }
    PD1.rowwise().normalize();
    PD2.rowwise().normalize();

    Eigen::MatrixXi isSingularity(N + 2, 1);
    Eigen::MatrixXi singularityIndex(N + 2, 1);
    igl::find_cross_field_singularities(V, F, PD1, PD2, isSingularity, singularityIndex, false);

    REQUIRE(isSingularity(0));
    REQUIRE(isSingularity(1));
    REQUIRE(singularityIndex(0) == 1);
    REQUIRE(singularityIndex(1) == 1);
    for (size_t i=2; i<N + 2; i++) {
        REQUIRE(!isSingularity(i));
        REQUIRE(singularityIndex(i) == 0);
    }
}
