
#include <test_common.h>
#include <igl/guess_extension.h>
#include <igl/pathinfo.h>
#include <cstdio>

class guess_extension : public ::testing::TestWithParam<std::string> {};

TEST_P(guess_extension, all_meshes)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::string path(test_common::data_path(GetParam()));
  // Load example mesh: GetParam() will be name of mesh file
  std::string d,b,e,f;
  igl::pathinfo(path,d,b,e,f);
  // Convert extension to lower case
  std::transform(e.begin(), e.end(), e.begin(), ::tolower);
  FILE * fp = fopen(path.c_str(),"r");
  std::string guess = igl::guess_extension(fp);
  ASSERT_EQ(guess,e);
}

INSTANTIATE_TEST_CASE_P
(
  all_meshes,
  guess_extension,
  ::testing::ValuesIn(test_common::all_meshes()),
  test_common::string_test_name
);
