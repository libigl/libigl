#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readMESH.h>
#include <fstream>

TEST_CASE("readMESH: single-tet","[igl]")
{
  std::ofstream("test.mesh")<< R"(MeshVersionFormatted 1
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
  igl::readMESH("test.mesh",V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}


TEST_CASE("readMESH: no-triangles-line","[igl]")
{
  std::ofstream("test.mesh")<< R"(MeshVersionFormatted 1
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
  igl::readMESH("test.mesh",V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}


TEST_CASE("readMESH: mesh-version-formatted-2","[igl]")
{
  std::ofstream("test.mesh")<< R"(MeshVersionFormatted 2
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
  igl::readMESH("test.mesh",V,T,F);
  REQUIRE(V.rows() == 4);
  REQUIRE(T.rows() == 1);
  REQUIRE(T(0,0) == 0);
  REQUIRE(F.rows() == 0);

}

