#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif

#include <test_common.h>
#include <igl/delaunay_triangulation.h>

namespace git_issue {

constexpr static double EPSILON_LENGTH = 0.0005; // Sketchup support 1/1000 precision
constexpr static double EPSILON_ANGLE = 0.0000000000005;
constexpr static double INV_EPSILON_LENGTH = 2000.0;

template<typename T>
inline int orient2d(const T &pa, const T &pb, const T &pc) {
    double acx, bcx, acy, bcy;
    acx = pa[0] - pc[0];
    bcx = pb[0] - pc[0];
    acy = pa[1] - pc[1];
    bcy = pb[1] - pc[1];

    double val = acx * bcy - acy * bcx;
    if (val < -EPSILON_LENGTH) {
        return -1;
    } else if (val > EPSILON_LENGTH) {
        return 1;
    } else {
        return 0;
    }
    return 0;
}

template <typename T>
inline int incircle(const T &pa, const T &pb, const T &pc, const T &pd) {
    double adx, ady, bdx, bdy, cdx, cdy;
    double abdet, bcdet, cadet;
    double alift, blift, clift;

    adx = pa[0] - pd[0];
    ady = pa[1] - pd[1];
    bdx = pb[0] - pd[0];
    bdy = pb[1] - pd[1];
    cdx = pc[0] - pd[0];
    cdy = pc[1] - pd[1];

    abdet = adx * bdy - bdx * ady;
    bcdet = bdx * cdy - cdx * bdy;
    cadet = cdx * ady - adx * cdy;
    alift = adx * adx + ady * ady;
    blift = bdx * bdx + bdy * bdy;
    clift = cdx * cdx + cdy * cdy;

    double val = alift * bcdet + blift * cadet + clift * abdet;
    if (val < -EPSILON_LENGTH) {
        return -1;
    } else if (val > EPSILON_LENGTH) {
        return 1;
    } else {
        return 0;
    }
    return 0;
}
}


TEST_CASE("delaunay_triangulation_issue_521", "[igl]") {
    using namespace Eigen;
    using namespace git_issue;
    MatrixXd V(16, 2);
    MatrixXi F;

    V << 4.55E-13, 2.33E-12,
        248.718, 249.939,
        463.602, 36.1059,
        764.953, 338.937,
        -1002.42, 2097.68,
        -1303.78, 1794.85,
        -1120.79, 1612.75,
        -1369.5, 1362.81,
        -1552.49, 1544.91,
        -1843.04, 1252.94,
        -75.6625, -505.806,
        214.883, -213.834,
        -191.721, 190.784,
        -1166.14, 1160.44,
        -917.424, 1410.38,
        56.9975, 440.723;

    const auto &orient2d_predicates = [](const double * pa, const double * pb, const double * pc) {
        return orient2d(pa, pb, pc);
    };
    const auto &incircle_predicates = [](const double * pa, const double * pb, const double * pc, const double * pd) {
        return incircle(pa, pb, pc, pd);
    };

    igl::delaunay_triangulation(V, orient2d_predicates, incircle_predicates, F);

    REQUIRE(F.rows() > 0);
    REQUIRE(F.cols() == 3);
    REQUIRE(F.maxCoeff() < 16);
    REQUIRE(F.minCoeff() == 0);
}

