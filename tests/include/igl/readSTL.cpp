#include <test_common.h>
#include <igl/readSTL.h>
#include <cstdio>
#include <fstream>

TEST_CASE("readSTL: ascii-long-word", "[igl]")
{
  const std::string filename = "readSTL_ascii-long-word.stl";
  std::ofstream(filename) <<
    "solid x\n" << std::string(250,'x') << "\nendsolid x\n";
  // readSTL closes fp
  FILE * fp = fopen(filename.c_str(),"rb");
  REQUIRE(fp != NULL);
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  REQUIRE(igl::readSTL(fp,V,F,N));
  REQUIRE(V.rows() == 0);
  REQUIRE(F.rows() == 0);
}
