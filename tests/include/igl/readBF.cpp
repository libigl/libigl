#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readBF.h>
#include <fstream>

TEST_CASE("readBF: valid-forest","[igl]")
{
  const std::string filename = "readBF_valid-forest.bf";
  // One root (weight -1, parent -1) followed by two non-root bones.
  std::ofstream(filename)<< R"(-1 -1 0 0 0
0 0 1 0 0
1 1 0 1 0
)";

  Eigen::VectorXi WI,bfP,P;
  Eigen::MatrixXd O,C;
  Eigen::MatrixXi BE;
  REQUIRE(igl::readBF(filename,WI,bfP,O,C,BE,P));
  REQUIRE(WI.rows() == 3);
  REQUIRE(C.rows() == 3);
  REQUIRE(BE.rows() == 2);
}

TEST_CASE("readBF: weight-index-out-of-range","[igl]")
{
  const std::string filename = "readBF_weight-index-out-of-range.bf";
  // A non-root bone whose weight index (5) exceeds the number of bones in the
  // forest. Without range checking this drives an out-of-bounds write into P
  // and BE.
  std::ofstream(filename)<< R"(-1 -1 0 0 0
5 0 1 1 1
)";

  Eigen::VectorXi WI,bfP,P;
  Eigen::MatrixXd O,C;
  Eigen::MatrixXi BE;
  REQUIRE_FALSE(igl::readBF(filename,WI,bfP,O,C,BE,P));
}

TEST_CASE("readBF: parent-index-out-of-range","[igl]")
{
  const std::string filename = "readBF_parent-index-out-of-range.bf";
  // A parent index (9) that references a joint that does not exist. Without
  // range checking this is read as an index into WI and dereferenced in
  // locate_tip.
  std::ofstream(filename)<< R"(-1 -1 0 0 0
0 9 1 1 1
)";

  Eigen::VectorXi WI,bfP,P;
  Eigen::MatrixXd O,C;
  Eigen::MatrixXi BE;
  REQUIRE_FALSE(igl::readBF(filename,WI,bfP,O,C,BE,P));
}

TEST_CASE("readBF: weighted-root","[igl]")
{
  const std::string filename = "readBF_weighted-root.bf";
  // A root (parent -1) with a weight index. Weighted bones look up their
  // parent's weight index, so without this check WI(-1) is read.
  std::ofstream(filename)<< R"(0 -1 0 0 0
-1 0 1 0 0
)";

  Eigen::VectorXi WI,bfP,P;
  Eigen::MatrixXd O,C;
  Eigen::MatrixXi BE;
  REQUIRE_FALSE(igl::readBF(filename,WI,bfP,O,C,BE,P));
}
