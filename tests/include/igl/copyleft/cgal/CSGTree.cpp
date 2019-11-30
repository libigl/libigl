#include <test_common.h>

#include <igl/copyleft/cgal/CSGTree.h>

TEST_CASE("CSGTree: extrusion", "[igl/copyleft/cgal]") {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path("extrusion.obj"), V, F);
    igl::copyleft::cgal::CSGTree tree(V, F);
    igl::copyleft::cgal::CSGTree inter(tree, tree, "i"); // returns error

    Eigen::MatrixXd V2 = inter.cast_V<Eigen::MatrixXd>();
    Eigen::MatrixXi F2 = inter.F();

    REQUIRE (V2.rows() == V.rows());
    REQUIRE (F2.rows() == F.rows());
}
