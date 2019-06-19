#include <test_common.h>
#include <igl/read_triangle_mesh.h>

TEST_CASE("read_triangle_mesh: cube_with_fold.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //test_common::load_mesh("cube_with_fold.ply", V, F);
    igl::read_triangle_mesh(test_common::data_path("cube_with_fold.ply"), V, F);
    REQUIRE (V.rows() == 26);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 48);
    REQUIRE (F.cols() == 3);
}

TEST_CASE("read_triangle_mesh: bunny.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path("bunny.ply"), V, F);
    REQUIRE (V.rows() == 34834);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 69451);
    REQUIRE (F.cols() == 3);
}
