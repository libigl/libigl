#include <test_common.h>
#include <igl/list_to_matrix.h>
#include <igl/STR.h>
#include <tuple>

namespace list_to_matrix
{
  typedef std::tuple<int/*n*/,int/*m*/> NM;
  // inline std::string NM_test_name(
  //   const ::testing::TestParamInfo<NM>& info)
  // {
  //   return STR(
  //     std::get<0>(info.param)<<"x"<<
  //     std::get<1>(info.param)<<"_");
  // };
}

TEST_CASE("ListToMatrixTest: matrix", "[igl]")
{
	const auto test_case = [](const list_to_matrix::NM &param)
	{
	  const int n = std::get<0>(param);
	  const int m = std::get<1>(param);
	  std::vector<std::vector<double> > vX(n,std::vector<double>(m));
	  for(int i = 0;i<n;i++)
	  {
	    for(int j = 0;j<m;j++)
	    {
	      vX[i][j] = i+j*n;
	    }
	  }
	  Eigen::MatrixXd mX;
	  igl::list_to_matrix(vX,mX);
	  for(int i = 0;i<n;i++)
	  {
	    for(int j = 0;j<m;j++)
	    {
	      REQUIRE (mX(i,j) == vX[i][j]);
	    }
	  }
	};

	std::vector<list_to_matrix::NM> params = {
	    list_to_matrix::NM{100,4},
	    list_to_matrix::NM{4,100},
	    list_to_matrix::NM{100,1},
	    list_to_matrix::NM{1,100},
    };

	test_common::run_test_cases(params, test_case);
}
