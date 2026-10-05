#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readMESH.h>
#include <cstdio>
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

TEST_CASE("readMESH: malformed-header","[igl]")
{
  std::string input;
  SECTION("empty file")
  {
    input = "";
  }
  SECTION("comments and whitespace only")
  {
    input = "# comment\n\n   \n";
  }
  SECTION("missing mesh version")
  {
    input = "MeshVersionFormatted\n";
  }
  SECTION("invalid mesh version")
  {
    input = "MeshVersionFormatted\ninvalid\n";
  }
  SECTION("missing dimension")
  {
    input = "MeshVersionFormatted 1\nDimension\n";
  }
  SECTION("invalid dimension")
  {
    input = "MeshVersionFormatted 1\nDimension\ninvalid\n";
  }

  const std::string filename = "readMESH_malformed-header.mesh";
  {
    std::ofstream file(filename);
    file << input;
    REQUIRE(file.good());
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  CHECK_FALSE(igl::readMESH(filename,V,T,F));
  CHECK(std::remove(filename.c_str()) == 0);
}

TEST_CASE("readMESH: split-header","[igl]")
{
  const std::string version = GENERATE("1","2");
  const std::string filename = "readMESH_split-header.mesh";
  {
    std::ofstream file(filename);
    file << "# comment\n\nMeshVersionFormatted\n" << version <<
      "\nDimension\n3\nVertices\n0\nTriangles\n0\nTetrahedra\n0\nEnd\n";
    REQUIRE(file.good());
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  CHECK(igl::readMESH(filename,V,T,F));
  CHECK(V.rows() == 0);
  CHECK(V.cols() == 3);
  CHECK(T.rows() == 0);
  CHECK(T.cols() == 4);
  CHECK(F.rows() == 0);
  CHECK(F.cols() == 3);
  CHECK(std::remove(filename.c_str()) == 0);
}
