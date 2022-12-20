#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readMESH.h>
#include <fstream>

TEST_CASE("readMESH: single-tet","[igl]")
{
  const std::string filename = "readMESH_single-tet.mesh";
  std::ofstream(filename)<< R"(MeshVersionFormatted 1
Dimension 3
Vertices
4
0 0 0 0
0 0 1 0
0 1 0 0
1 0 0 0
Triangles
0
Tetrahedra
1
1 2 3 4 0
End
  )";

  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  igl::readMESH(filename,V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}


TEST_CASE("readMESH: no-triangles-line","[igl]")
{
  const std::string filename = "readMESH_no-triangles-line.mesh";
  std::ofstream(filename)<< R"(MeshVersionFormatted 1
Dimension 3
Vertices
4
0 0 0 0
0 0 1 0
0 1 0 0
1 0 0 0
Tetrahedra
1
1 2 3 4 0
  )";

  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  igl::readMESH(filename,V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}


TEST_CASE("readMESH: mesh-version-formatted-2","[igl]")
{
  const std::string filename = "readMESH_mesh-version-formatted-2.mesh";
  std::ofstream(filename)<< R"(MeshVersionFormatted 2
Dimension 3
Vertices
4
0 0 0 0
0 0 1 0
0 1 0 0
1 0 0 0
Triangles
0
Tetrahedra
1
1 2 3 4 0
End
  )";

  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  igl::readMESH(filename,V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}

