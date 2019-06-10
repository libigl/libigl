#include <test_common.h>
#include <igl/sort.h>
#include <igl/STR.h>
#include <tuple>

namespace sort
{
  typedef std::tuple<int/*n*/,int/*m*/,int/*dim*/,bool/*ascending*/> NMDimAscending;
  // inline std::string NMDimAscending_test_name(
  //   const ::testing::TestParamInfo<NMDimAscending>& info)
  // {
  //   return STR(
  //     std::get<0>(info.param)<<"x"<<
  //     std::get<1>(info.param)<<"_"<<
  //     "dim_"<<std::get<2>(info.param)<<"_"<<
  //     "ascending_"<<(std::get<3>(info.param)?"true":"false"));
  // };
}

TEST_CASE("SortTest: random", "[igl]")
{
	const auto test_case = [](const sort::NMDimAscending &param)
	{
		const int n = std::get<0>(param);
		const int m = std::get<1>(param);
		const int dim = std::get<2>(param);
		const bool ascending = std::get<3>(param);
		Eigen::MatrixXd X = Eigen::MatrixXd::Random(n,m);
  		// sort ascending
		Eigen::MatrixXd Y;
		Eigen::MatrixXi IX;
		igl::sort(X,dim,ascending,Y,IX);
		REQUIRE (Y.rows() == X.rows());
		REQUIRE (Y.cols() == X.cols());
		REQUIRE (IX.rows() == X.rows());
		REQUIRE (IX.cols() == X.cols());
		for(int i = 0;i<n;i++)
		{
			for(int j = 0;j<m;j++)
			{
				REQUIRE (X(dim==1?IX(i,j):i,dim==2?IX(i,j):j) == Y(i,j));
			}
		}
		for(int i = (dim==1?1:0);i<n;i++)
		{
			for(int j = (dim==2?1:0);j<m;j++)
			{
				if(ascending)
				{
					REQUIRE (Y(i,j) >= Y(i-(dim==1?1:0),j-(dim==2?1:0)));
				}else
				{
					REQUIRE (Y(i,j) <= Y(i-(dim==1?1:0),j-(dim==2?1:0)));
				}
			}
		}
	};

	std::vector<sort::NMDimAscending> params = {
		sort::NMDimAscending{100,3,1,true},
		sort::NMDimAscending{100,3,2,true},
		sort::NMDimAscending{100,3,1,false},
		sort::NMDimAscending{100,3,2,false},
		sort::NMDimAscending{3,100,1,true},
		sort::NMDimAscending{3,100,2,true},
		sort::NMDimAscending{3,100,1,false},
		sort::NMDimAscending{3,100,2,false},
		sort::NMDimAscending{100,2,1,true},
		sort::NMDimAscending{100,2,2,true},
		sort::NMDimAscending{100,2,1,false},
		sort::NMDimAscending{100,2,2,false},
		sort::NMDimAscending{2,100,1,true},
		sort::NMDimAscending{2,100,2,true},
		sort::NMDimAscending{2,100,1,false},
		sort::NMDimAscending{2,100,2,false},
		sort::NMDimAscending{100,4,1,true},
		sort::NMDimAscending{100,4,2,true},
		sort::NMDimAscending{100,4,1,false},
		sort::NMDimAscending{100,4,2,false},
		sort::NMDimAscending{4,100,1,true},
		sort::NMDimAscending{4,100,2,true},
		sort::NMDimAscending{4,100,1,false},
		sort::NMDimAscending{4,100,2,false},
	};

	test_common::run_test_cases(params, test_case);
}
