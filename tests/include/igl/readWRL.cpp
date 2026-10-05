#include <test_common.h>
#include <igl/readWRL.h>
#include <cstdio>
#include <fstream>

TEST_CASE("readWRL: missing-coordIndex", "[igl]")
{
  const std::string filename = "readWRL_missing-coordIndex.wrl";
  {
    std::ofstream file(filename);
    file << R"(#VRML V2.0 utf8
Shape {
  geometry IndexedFaceSet {
    coord Coordinate {
      point [
        0 0 0,
        1 0 0,
        0 1 0
      ]
    }
  }
}
)";
    REQUIRE(file.good());
  }
  std::vector<std::vector<double> > V;
  std::vector<std::vector<int> > F;
  CHECK_FALSE(igl::readWRL(filename,V,F));
  CHECK(std::remove(filename.c_str()) == 0);
}

TEST_CASE("readWRL: simple", "[igl]")
{
  const std::string filename = "readWRL_simple.wrl";
  {
    std::ofstream file(filename);
    file << R"(#VRML V2.0 utf8
Shape {
  geometry IndexedFaceSet {
    coord Coordinate {
      point [
        0 0 0,
        1 0 0,
        0 1 0
      ]
    }
    coordIndex [
      0,1,2,-1
    ]
  }
}
)";
    REQUIRE(file.good());
  }
  std::vector<std::vector<double> > V;
  std::vector<std::vector<int> > F;
  CHECK(igl::readWRL(filename,V,F));
  const std::vector<std::vector<double> > expected_V =
    {{0,0,0},{1,0,0},{0,1,0}};
  const std::vector<std::vector<int> > expected_F = {{0,1,2}};
  CHECK(V == expected_V);
  CHECK(F == expected_F);
  CHECK(std::remove(filename.c_str()) == 0);
}
