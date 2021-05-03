#include <test_common.h>
#include <igl/predicates/predicates.h>
#include <igl/triangle/triangulate.h>
#include <limits>

TEST_CASE("predicates and triangle", "[igl][predicates][triangle]") {
    using namespace igl::predicates;
    using Scalar = double;
    igl::predicates::exactinit();

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
