#include <test_common.h>
#include <igl/readOFF.h>
#include <fstream>

TEST_CASE("readOFF: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOFF
    // directly...
    igl::read_triangle_mesh(test_common::data_path("cube.off"), V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 12);
    REQUIRE (F.cols() == 3);
}

TEST_CASE("readOFF: long-header", "[igl]")
{
  const std::string filename = "readOFF_long-header.off";
  std::ofstream(filename) << std::string(2000,'x') <<
    "\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n";
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  REQUIRE_FALSE(igl::readOFF(filename,V,F));
}

TEST_CASE("readOFF: long-comments", "[igl]")
{
  const std::string filename = "readOFF_long-comments.off";
  const std::string comment = "#" + std::string(2000,'x') + "\n";
  std::ofstream(filename) <<
    "OFF\n3 1 0\n0 0 0\n# a\n" << comment << "1 0 0\n0 1 0\n" <<
    comment << "3 0 1 2\n";
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  REQUIRE(igl::readOFF(filename,V,F));
  REQUIRE(V.rows() == 3);
  REQUIRE(F.rows() == 1);
  REQUIRE(F(0,2) == 2);
}
