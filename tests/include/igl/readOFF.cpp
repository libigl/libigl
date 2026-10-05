#include <test_common.h>
#include <igl/readOFF.h>
#include <cstdio>
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

TEST_CASE("readOFF: malformed-input", "[igl]")
{
  const std::string vertices = "OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n";
  std::string input;
  SECTION("empty file")
  {
    input = "";
  }
  SECTION("missing counts")
  {
    input = "OFF\n";
  }
  SECTION("comments without counts")
  {
    input = "OFF\n# comment\n";
  }
  SECTION("incomplete counts")
  {
    input = "OFF\n3 1\n";
  }
  SECTION("invalid counts")
  {
    input = "OFF\ninvalid 1 0\n";
  }
  SECTION("missing vertices")
  {
    input = "OFF\n1 0 0\n";
  }
  SECTION("truncated vertices")
  {
    input = "OFF\n2 0 0\n0 0 0\n";
  }
  SECTION("missing faces")
  {
    input = vertices;
  }
  SECTION("truncated intermediate face index")
  {
    input = vertices + "3 0";
  }
  SECTION("truncated last face index")
  {
    input = vertices + "3 0 1\n";
  }
  SECTION("invalid intermediate face index")
  {
    input = vertices + "3 0 invalid 2\n";
  }
  SECTION("invalid last face index")
  {
    input = vertices + "3 0 1 invalid\n";
  }
  SECTION("unterminated comment without a face")
  {
    input = vertices + "# comment";
  }

  const std::string filename = "readOFF_malformed-input.off";
  {
    std::ofstream file(filename);
    file << input;
    REQUIRE(file.good());
  }
  std::vector<std::vector<double> > V,N,C;
  std::vector<std::vector<int> > F;
  CHECK_FALSE(igl::readOFF(filename,V,F,N,C));
  CHECK(std::remove(filename.c_str()) == 0);
}

TEST_CASE("readOFF: comments", "[igl]")
{
  std::string comment;
  SECTION("empty comment")
  {
    comment = "#\n";
  }
  SECTION("nonempty comment")
  {
    comment = "# comment\n";
  }
  SECTION("multiple comment markers")
  {
    comment = "### comment\n";
  }
  SECTION("long comment")
  {
    comment = "#" + std::string(3000,'x') + "\n";
  }

  const std::string filename = "readOFF_comments.off";
  {
    std::ofstream file(filename);
    file << "OFF\n# header comment\n3 2 0\n#\n" << comment <<
      "0 0 0\n1 0 0\n0 1 0\n" << comment <<
      "3 0 1 2\n" << comment << "3 2 1 0\n";
    REQUIRE(file.good());
  }
  std::vector<std::vector<double> > V,N,C;
  std::vector<std::vector<int> > F;
  CHECK(igl::readOFF(filename,V,F,N,C));
  const std::vector<std::vector<double> > expected_V =
    {{0,0,0},{1,0,0},{0,1,0}};
  const std::vector<std::vector<int> > expected_F = {{0,1,2},{2,1,0}};
  CHECK(V == expected_V);
  CHECK(F == expected_F);
  CHECK(N.empty());
  CHECK(C.empty());
  CHECK(std::remove(filename.c_str()) == 0);
}
