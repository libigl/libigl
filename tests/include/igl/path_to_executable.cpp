#include <test_common.h>
#include <igl/path_to_executable.h>

#include <iostream>


TEST_CASE("path_to_executable: example", "[igl]")
{
  std::string path_to_executable = igl::path_to_executable();
  REQUIRE(0 < path_to_executable.size());
  // check if path_to_executable ends with correct file name
  std::string ending = "libigl_tests";
  REQUIRE(ending.size() < path_to_executable.size());
  REQUIRE(std::equal(ending.rbegin(), ending.rend(), path_to_executable.rbegin()));
}
