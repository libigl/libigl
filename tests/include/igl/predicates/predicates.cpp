#include <test_common.h>
#include <igl/predicates/orient2d.h>
#include <igl/predicates/orient3d.h>
#include <igl/predicates/incircle.h>
#include <igl/predicates/insphere.h>
#include <igl/predicates/exactinit.h>
#include <limits>

// Didn't have the stamina to break the tests into separate files but they
// should be
TEST_CASE("predicates", "[igl][predicates]") {
    using namespace igl::predicates;
    using Scalar = double;
    igl::predicates::exactinit();

    SECTION("2D") {
        using Point = Eigen::Matrix<Scalar, 2, 1>;
        Point a(2,1),b(2,1),c(2,1),d(2,1),e(2,1),f(2,1);
        a << 0.0, 0.0;
        b << 1.0, 0.0;
        c << 1.0, 1.0;
        d << 2.0, 0.0;
        e << 0.0, 1.0;
        f << 0.5, 0.5;

        REQUIRE(orient2d(a, b, c) == igl::Orientation::POSITIVE);
        REQUIRE(orient2d(a, c, b) == igl::Orientation::NEGATIVE);
        REQUIRE(orient2d(a, b, b) == igl::Orientation::COLLINEAR);
        REQUIRE(orient2d(a, a, a) == igl::Orientation::COLLINEAR);
        REQUIRE(orient2d(a, b, d) == igl::Orientation::COLLINEAR);
        REQUIRE(orient2d(a, f, c) == igl::Orientation::COLLINEAR);

        REQUIRE(incircle(a,b,c,e) == igl::Orientation::COCIRCULAR);
        REQUIRE(incircle(a,b,c,a) == igl::Orientation::COCIRCULAR);
        REQUIRE(incircle(a,b,c,d) == igl::Orientation::OUTSIDE);
        REQUIRE(incircle(a,b,c,f) == igl::Orientation::INSIDE);
    }

    SECTION("3D") {
        using Point = Eigen::Matrix<Scalar, 3, 1>;
        Point a(3,1),b(3,1),c(3,1),d(3,1),e(3,1),f(3,1);
        a << 0.0, 0.0, 0.0;
        b << 1.0, 0.0, 0.0;
        c << 0.0, 1.0, 0.0;
        d << 0.0, 0.0, 1.0;
        e << 1.0, 1.0, 1.0;
        f << std::numeric_limits<Scalar>::epsilon(), 0.0, 0.0;

        REQUIRE(orient3d(a, b, c, d) == igl::Orientation::NEGATIVE);
        REQUIRE(orient3d(a, b, d, c) == igl::Orientation::POSITIVE);
        REQUIRE(orient3d(a, b, d, d) == igl::Orientation::COPLANAR);
        REQUIRE(orient3d(a, a, a, a) == igl::Orientation::COPLANAR);
        REQUIRE(orient3d(a, b, f, c) == igl::Orientation::COPLANAR);

        REQUIRE(insphere(a, b, c, d, e)                  == igl::Orientation::COSPHERICAL);
        REQUIRE(insphere(a, b, d, e, c)                  == igl::Orientation::COSPHERICAL);
        REQUIRE(insphere(b, c, e, d, ((a+b)*0.5).eval()) == igl::Orientation::INSIDE);
        REQUIRE(insphere(b, c, e, d, (-f).eval())        == igl::Orientation::OUTSIDE);
        REQUIRE(insphere(f, b, d, c, e)                  == igl::Orientation::INSIDE);
    }
}
