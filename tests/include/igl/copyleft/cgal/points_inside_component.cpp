#include <test_common.h>

#include <igl/copyleft/cgal/points_inside_component.h>
#include <igl/upsample.h>
#include <limits>

TEST_CASE("PointInsideComponent: simple", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V1, F1);

    Eigen::MatrixXd P(4, 3);
    P << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    Eigen::VectorXi inside;

    CHECK_NOTHROW (igl::copyleft::cgal::points_inside_component(V1, F1, P, inside));
    REQUIRE (inside[0] == 1);
    REQUIRE (inside[1] == 0);
    REQUIRE (inside[2] == 0);
    REQUIRE (inside[3] == 0);
}

TEST_CASE("PointInsideComponent: near_boundary", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V1, F1);

    const double EPS = std::numeric_limits<double>::epsilon();
    Eigen::MatrixXd P(6, 3);
    P << 0.5 + EPS, 0.0, 0.0,
         0.0, 0.5 + EPS, 0.0,
         0.0, 0.0, 0.5 + EPS,
         0.5 - EPS, 0.0, 0.0,
         0.0, 0.5 - EPS, 0.0,
         0.0, 0.0, 0.5 - EPS;

    Eigen::VectorXi inside;
    CHECK_NOTHROW (igl::copyleft::cgal::points_inside_component(V1, F1, P, inside));
    REQUIRE (inside[0] == 0);
    REQUIRE (inside[1] == 0);
    REQUIRE (inside[2] == 0);
    REQUIRE (inside[3] == 1);
    REQUIRE (inside[4] == 1);
    REQUIRE (inside[5] == 1);
}

TEST_CASE("PointInsideComponent: near_corner", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V1, F1);

    const double EPS = std::numeric_limits<double>::epsilon();
    Eigen::MatrixXd P_out(8, 3);
    P_out << 0.5 + EPS, 0.5 + EPS, 0.5 + EPS,
            -0.5 - EPS, 0.5 + EPS, 0.5 + EPS,
             0.5 + EPS,-0.5 - EPS, 0.5 + EPS,
            -0.5 - EPS,-0.5 - EPS, 0.5 + EPS,
             0.5 + EPS, 0.5 + EPS,-0.5 - EPS,
            -0.5 - EPS, 0.5 + EPS,-0.5 - EPS,
             0.5 + EPS,-0.5 - EPS,-0.5 - EPS,
            -0.5 - EPS,-0.5 - EPS,-0.5 - EPS;

    Eigen::VectorXi inside;
    CHECK_NOTHROW (igl::copyleft::cgal::points_inside_component(V1, F1, P_out, inside));
    REQUIRE ((inside.array()==0).all());

    Eigen::MatrixXd P_in(8, 3);
    P_in << 0.5 - EPS, 0.5 - EPS, 0.5 - EPS,
           -0.5 + EPS, 0.5 - EPS, 0.5 - EPS,
            0.5 - EPS,-0.5 + EPS, 0.5 - EPS,
           -0.5 + EPS,-0.5 + EPS, 0.5 - EPS,
            0.5 - EPS, 0.5 - EPS,-0.5 + EPS,
           -0.5 + EPS, 0.5 - EPS,-0.5 + EPS,
            0.5 - EPS,-0.5 + EPS,-0.5 + EPS,
           -0.5 + EPS,-0.5 + EPS,-0.5 + EPS;
    CHECK_NOTHROW (igl::copyleft::cgal::points_inside_component(V1, F1, P_in, inside));
    REQUIRE ((inside.array()==1).all());
}

TEST_CASE("PointInsideComponent: closest_vertex_with_collinear_neighbors", "[igl/copyleft/cgal]")
{
    // Upsampling the cube puts a vertex at the midpoint of each cube edge.
    // The one-ring of such a vertex contains the two endpoints of that cube
    // edge, which are collinear with it. Queries whose closest point is such a
    // midpoint vertex must not trip over the degenerate (collinear) candidate
    // separating plane.
    Eigen::MatrixXd V0, V;
    Eigen::MatrixXi F0, F;
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V0, F0);
    igl::upsample(V0, F0, V, F);

    // Push each cube-edge midpoint outward along its bisector; the closest
    // point on the surface is exactly that (mesh) vertex.
    Eigen::MatrixXd P(12, 3);
    int k = 0;
    for (int a = 0; a < 3; a++)
    {
        const int b = (a+1)%3;
        for (int sa : {-1, 1}) for (int sb : {-1, 1})
        {
            P.row(k).setZero();
            P(k, a) = sa * 1.0;
            P(k, b) = sb * 1.0;
            k++;
        }
    }
    Eigen::VectorXi inside;
    REQUIRE_NOTHROW (igl::copyleft::cgal::points_inside_component(V, F, P, inside));
    REQUIRE ((inside.array()==0).all());
}

TEST_CASE("PointInsideComponent: inconsistent_edge_orientation_throws", "[igl/copyleft/cgal]")
{
    // Flip one side of the cube so the edges bounding it are traversed in the
    // same direction by both incident faces. A query whose closest point lies
    // on such an edge cannot be classified and must be rejected rather than
    // silently decided by the first face alone.
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
    // Flip every face lying in the plane x = 0.5.
    for (int f = 0; f < F.rows(); f++)
    {
        bool on_side = true;
        for (int c = 0; c < 3; c++) on_side &= V(F(f,c),0) == 0.5;
        if (on_side) std::swap(F(f,1), F(f,2));
    }
    Eigen::MatrixXd P(1, 3);
    P << 1.0, 1.0, 0.0; // closest to the interior of the edge x=y=0.5
    Eigen::VectorXi inside;
    REQUIRE_THROWS (igl::copyleft::cgal::points_inside_component(V, F, P, inside));
}
