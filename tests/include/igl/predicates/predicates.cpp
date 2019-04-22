#include <test_common.h>
#include <igl/predicates/predicates.h>
#include <limits>
#include <igl/triangle/triangulate.h>

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

        REQUIRE(orient2d(a, b, c) == Orientation::POSITIVE);
        REQUIRE(orient2d(a, c, b) == Orientation::NEGATIVE);
        REQUIRE(orient2d(a, b, b) == Orientation::COLLINEAR);
        REQUIRE(orient2d(a, a, a) == Orientation::COLLINEAR);
        REQUIRE(orient2d(a, b, d) == Orientation::COLLINEAR);
        REQUIRE(orient2d(a, f, c) == Orientation::COLLINEAR);

        REQUIRE(incircle(a,b,c,e) == Orientation::COCIRCULAR);
        REQUIRE(incircle(a,b,c,a) == Orientation::COCIRCULAR);
        REQUIRE(incircle(a,b,c,d) == Orientation::OUTSIDE);
        REQUIRE(incircle(a,b,c,f) == Orientation::INSIDE);
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

        REQUIRE(orient3d(a, b, c, d) == Orientation::NEGATIVE);
        REQUIRE(orient3d(a, b, d, c) == Orientation::POSITIVE);
        REQUIRE(orient3d(a, b, d, d) == Orientation::COPLANAR);
        REQUIRE(orient3d(a, a, a, a) == Orientation::COPLANAR);
        REQUIRE(orient3d(a, b, f, c) == Orientation::COPLANAR);

        REQUIRE(insphere(a, b, c, d, e) == Orientation::COSPHERICAL);
        REQUIRE(insphere(a, b, d, e, c) == Orientation::COSPHERICAL);
        REQUIRE(insphere(b, c, e, d, ((a+b)*0.5).eval()) == Orientation::INSIDE);
        REQUIRE(insphere(b, c, e, d, (-f).eval()) == Orientation::OUTSIDE);
        REQUIRE(insphere(f, b, d, c, e) == Orientation::INSIDE);
    }

    SECTION("Predicate and triangle") {
        Eigen::Matrix<double, -1, -1> vertices(4, 2);
        Eigen::Matrix<double, -1, -1> holes;
        Eigen::Matrix<int, -1, -1> edges;
        vertices << 0.0, 0.0,
                    1.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0;

        Eigen::Matrix<double, -1, -1> out_vertices;
        Eigen::Matrix<int, -1, -1> out_faces;

        // Run constrained Delaunay.
        igl::triangle::triangulate(vertices, edges, holes, "QcYY",
                out_vertices, out_faces);
        REQUIRE(out_vertices.rows() == 4);
        REQUIRE(out_vertices.cols() == 2);
        REQUIRE(out_faces.rows() == 2);
        REQUIRE(out_faces.cols() == 3);
    }
}
