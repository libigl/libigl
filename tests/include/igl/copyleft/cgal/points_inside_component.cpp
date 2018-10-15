#include <test_common.h>

#include <igl/copyleft/cgal/points_inside_component.h>
#include <limits>

namespace PointInsideComponentHelper {

TEST(PointInsideComponent, simple) {
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    test_common::load_mesh("cube.obj", V1, F1);

    Eigen::MatrixXd P(4, 3);
    P << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    Eigen::VectorXi inside;

    EXPECT_NO_THROW(igl::copyleft::cgal::points_inside_component(V1, F1, P, inside));
    ASSERT_EQ(1, inside[0]);
    ASSERT_EQ(0, inside[1]);
    ASSERT_EQ(0, inside[2]);
    ASSERT_EQ(0, inside[3]);
}

TEST(PointInsideComponent, near_boundary) {
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    test_common::load_mesh("cube.obj", V1, F1);

    const double EPS = std::numeric_limits<double>::epsilon();
    Eigen::MatrixXd P(6, 3);
    P << 0.5 + EPS, 0.0, 0.0,
         0.0, 0.5 + EPS, 0.0,
         0.0, 0.0, 0.5 + EPS,
         0.5 - EPS, 0.0, 0.0,
         0.0, 0.5 - EPS, 0.0,
         0.0, 0.0, 0.5 - EPS;

    Eigen::VectorXi inside;
    EXPECT_NO_THROW(igl::copyleft::cgal::points_inside_component(V1, F1, P, inside));
    ASSERT_EQ(0, inside[0]);
    ASSERT_EQ(0, inside[1]);
    ASSERT_EQ(0, inside[2]);
    ASSERT_EQ(1, inside[3]);
    ASSERT_EQ(1, inside[4]);
    ASSERT_EQ(1, inside[5]);
}

TEST(PointInsideComponent, near_corner) {
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    test_common::load_mesh("cube.obj", V1, F1);

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
    EXPECT_NO_THROW(igl::copyleft::cgal::points_inside_component(V1, F1, P_out, inside));
    ASSERT_TRUE((inside.array()==0).all());

    Eigen::MatrixXd P_in(8, 3);
    P_in << 0.5 - EPS, 0.5 - EPS, 0.5 - EPS,
           -0.5 + EPS, 0.5 - EPS, 0.5 - EPS,
            0.5 - EPS,-0.5 + EPS, 0.5 - EPS,
           -0.5 + EPS,-0.5 + EPS, 0.5 - EPS,
            0.5 - EPS, 0.5 - EPS,-0.5 + EPS,
           -0.5 + EPS, 0.5 - EPS,-0.5 + EPS,
            0.5 - EPS,-0.5 + EPS,-0.5 + EPS,
           -0.5 + EPS,-0.5 + EPS,-0.5 + EPS;
    EXPECT_NO_THROW(igl::copyleft::cgal::points_inside_component(V1, F1, P_in, inside));
    ASSERT_TRUE((inside.array()==1).all());
}

}
