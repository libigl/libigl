#include <test_common.h>
#include <igl/readPLY.h>
#include <fstream>
#include <string>
#include <vector>

TEST_CASE("readPLY: cube_with_fold.ply", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    REQUIRE (igl::readPLY(test_common::data_path("cube_with_fold.ply"), V, F));
    REQUIRE (V.rows() == 26);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 48);
    REQUIRE (F.cols() == 3);
}

TEST_CASE("readPLY: bunny.ply", "[igl]")
{
    std::ifstream f(test_common::data_path("bunny.ply"));
    REQUIRE (f.good());
    f.close();
    
    Eigen::MatrixXd V,N,UV,VD,FD,ED;
    std::vector<std::string> Vheader,Fheader,Eheader,comments;

    Eigen::MatrixXi F,E;

    REQUIRE (igl::readPLY(test_common::data_path("bunny.ply"), V, F, E, N, UV, VD,Vheader, FD,Fheader, ED,Eheader,comments));
    REQUIRE (V.rows() == 35947);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 69451);
    REQUIRE (F.cols() == 3);
    // no edge data
    REQUIRE (E.rows() == 0);
    REQUIRE (E.cols() == 0);

    // no normals or texture coordinates
    REQUIRE (N.rows() == 0);
    REQUIRE (N.cols() == 0);
    REQUIRE (UV.rows() == 0);
    REQUIRE (UV.cols() == 0);

    // this bunny have additonal data
    REQUIRE (VD.rows() == 35947);
    REQUIRE (VD.cols() == 2);

    REQUIRE (Vheader.size() == 2);
    REQUIRE (Vheader[0] == "confidence" );
    REQUIRE (Vheader[1] == "intensity" );

    // no Face data or edge data
    REQUIRE (FD.rows() == 0);
    REQUIRE (FD.cols() == 0);
    REQUIRE (Fheader.size() == 0);

    REQUIRE (ED.rows() == 0);
    REQUIRE (ED.cols() == 0);
    REQUIRE (Eheader.size() == 0);

    // there are comments
    REQUIRE (comments.size() == 2);
}

TEST_CASE("readPLY: mesh_error.ply", "[igl]")
{   
    // test on a non-existent file
    std::ifstream f(test_common::data_path("mesh_error.ply"));
    REQUIRE (f.good() == false);
    f.close();

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    REQUIRE (igl::readPLY(test_common::data_path("mesh_error.ply"), V, F) == false);
    REQUIRE (V.rows() == 0);
    REQUIRE (F.rows() == 0);
}

TEST_CASE("readPLY: quad_cube.ply", "[igl]")
{
    // small qube from blender
    const char *ply_quad_cube=
"ply\n"
"format ascii 1.0\n"
"comment Created by Blender 2.81 (sub 16) - www.blender.org, source file: ''\n"
"element vertex 24\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property float nx\n"
"property float ny\n"
"property float nz\n"
"property float s\n"
"property float t\n"
"element face 6\n"
"property list uchar uint vertex_indices\n"
"end_header\n"
"1.000000 1.000000 1.000000 0.000000 -0.000000 1.000000 0.625000 0.500000\n"
"-1.000000 1.000000 1.000000 0.000000 -0.000000 1.000000 0.875000 0.500000\n"
"-1.000000 -1.000000 1.000000 0.000000 -0.000000 1.000000 0.875000 0.750000\n"
"1.000000 -1.000000 1.000000 0.000000 -0.000000 1.000000 0.625000 0.750000\n"
"1.000000 -1.000000 -1.000000 0.000000 -1.000000 0.000000 0.375000 0.750000\n"
"1.000000 -1.000000 1.000000 0.000000 -1.000000 0.000000 0.625000 0.750000\n"
"-1.000000 -1.000000 1.000000 0.000000 -1.000000 0.000000 0.625000 1.000000\n"
"-1.000000 -1.000000 -1.000000 0.000000 -1.000000 0.000000 0.375000 1.000000\n"
"-1.000000 -1.000000 -1.000000 -1.000000 -0.000000 0.000000 0.375000 0.000000\n"
"-1.000000 -1.000000 1.000000 -1.000000 -0.000000 0.000000 0.625000 0.000000\n"
"-1.000000 1.000000 1.000000 -1.000000 -0.000000 0.000000 0.625000 0.250000\n"
"-1.000000 1.000000 -1.000000 -1.000000 -0.000000 0.000000 0.375000 0.250000\n"
"-1.000000 1.000000 -1.000000 0.000000 0.000000 -1.000000 0.125000 0.500000\n"
"1.000000 1.000000 -1.000000 0.000000 0.000000 -1.000000 0.375000 0.500000\n"
"1.000000 -1.000000 -1.000000 0.000000 0.000000 -1.000000 0.375000 0.750000\n"
"-1.000000 -1.000000 -1.000000 0.000000 0.000000 -1.000000 0.125000 0.750000\n"
"1.000000 1.000000 -1.000000 1.000000 -0.000000 0.000000 0.375000 0.500000\n"
"1.000000 1.000000 1.000000 1.000000 -0.000000 0.000000 0.625000 0.500000\n"
"1.000000 -1.000000 1.000000 1.000000 -0.000000 0.000000 0.625000 0.750000\n"
"1.000000 -1.000000 -1.000000 1.000000 -0.000000 0.000000 0.375000 0.750000\n"
"-1.000000 1.000000 -1.000000 0.000000 1.000000 0.000000 0.375000 0.250000\n"
"-1.000000 1.000000 1.000000 0.000000 1.000000 0.000000 0.625000 0.250000\n"
"1.000000 1.000000 1.000000 0.000000 1.000000 0.000000 0.625000 0.500000\n"
"1.000000 1.000000 -1.000000 0.000000 1.000000 0.000000 0.375000 0.500000\n"
"4 0 1 2 3\n"
"4 4 5 6 7\n"
"4 8 9 10 11\n"
"4 12 13 14 15\n"
"4 16 17 18 19\n"
"4 20 21 22 23\n";

    std::ofstream f("quad_cube.ply");
    f.write(ply_quad_cube,strlen(ply_quad_cube));
    f.close();

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    REQUIRE (igl::readPLY("quad_cube.ply", V, F));
    REQUIRE (V.rows() == 24);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 6);
    REQUIRE (F.cols() == 4);
}