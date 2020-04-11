#include <test_common.h>

#include <igl/copyleft/cgal/points_inside_component.h>
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
