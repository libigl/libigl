#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readNODE.h>
#include <fstream>

TEST_CASE("readNODE: single-point","[igl]")
{
  const std::string filename = "readNODE_single-point.node";
  std::ofstream(filename)<< R"(1 3 0 0
1 0.5 0.25 0.75
)";

  Eigen::MatrixXd V;
  Eigen::MatrixXi I;
  REQUIRE(igl::readNODE(filename,V,I));
  REQUIRE(V.rows() == 1);
  REQUIRE(V.cols() == 3);
  REQUIRE(I.rows() == 1);
  REQUIRE(V(0,0) == Approx(0.5));
  REQUIRE(V(0,1) == Approx(0.25));
  REQUIRE(V(0,2) == Approx(0.75));
}

TEST_CASE("readNODE: more-points-than-declared","[igl]")
{
  // The header declares a single point but the body lists several. The reader
  // must reject this rather than writing V(p,.)/I(p) past the allocated rows.
  const std::string filename = "readNODE_more-points-than-declared.node";
  std::ofstream(filename)<< R"(1 3 0 0
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
)";

  Eigen::MatrixXd V;
  Eigen::MatrixXi I;
  REQUIRE_FALSE(igl::readNODE(filename,V,I));
}
