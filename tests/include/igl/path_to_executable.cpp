#include <test_common.h>
#include <igl/path_to_executable.h>

#include <iostream>


TEST_CASE("path_to_executable: example", "[igl]")
{
  std::string path_to_executable = igl::path_to_executable();
  REQUIRE(0 < path_to_executable.size());
  // check if path_to_executable ends with correct file name, on windows .exe suffix is added.
  std::string executable = "libigl_tests";
  int pos = path_to_executable.length()-(executable.length() + 4/*".exe"*/);
  REQUIRE( std::string::npos != path_to_executable.find(executable, pos));
}
