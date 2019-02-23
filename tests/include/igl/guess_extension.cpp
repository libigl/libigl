
#include <test_common.h>
#include <igl/guess_extension.h>
#include <igl/pathinfo.h>
#include <cstdio>


TEST_CASE("guess_extension: all_meshes", "[igl]")
{
	const auto test_case = [](const std::string &param)
	{
	  Eigen::MatrixXd V;
	  Eigen::MatrixXi F;
	  std::string path(test_common::data_path(param));
	  // Load example mesh: GetParam() will be name of mesh file
	  std::string d,b,e,f;
	  igl::pathinfo(path,d,b,e,f);
	  // Convert extension to lower case
	  std::transform(e.begin(), e.end(), e.begin(), ::tolower);
	  FILE * fp = fopen(path.c_str(),"r");
	  std::string guess = igl::guess_extension(fp);
	  REQUIRE (e == guess);
	};

	test_common::run_test_cases(test_common::all_meshes(), test_case);
}
