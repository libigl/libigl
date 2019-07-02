#include <test_common.h>
#include <igl/readPLY.h>

TEST_CASE("readPLY: cube_with_fold.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //test_common::load_mesh("cube_with_fold.ply", V, F);
    REQUIRE (igl::readPLY(test_common::data_path("cube_with_fold.ply"), V, F));
    REQUIRE (V.rows() == 26);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 48);
    REQUIRE (F.cols() == 3);
}

TEST_CASE("readPLY: bunny.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    REQUIRE (igl::readPLY(test_common::data_path("bunny.ply"), V, F));
    REQUIRE (V.rows() == 34834);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 69451);
    REQUIRE (F.cols() == 3);
}

TEST_CASE("readPLY: mesh_error.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    REQUIRE (igl::readPLY(test_common::data_path("mesh_error.ply"), V, F) == false);
    REQUIRE (V.rows() == 0);
    REQUIRE (F.rows() == 0);
}
