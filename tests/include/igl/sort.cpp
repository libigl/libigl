#include <test_common.h>
#include <igl/sort.h>
#include <igl/sortrows.h>
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

TEST_CASE("SortTest: stability", "[igl]")
{
	Eigen::VectorXi X(6);
	X << 2, 1, 2, 1, 2, 1;
	Eigen::VectorXi Y, IX;

	Eigen::VectorXi expected_asc_ix(6);
	expected_asc_ix << 1, 3, 5, 0, 2, 4;
	Eigen::VectorXi expected_desc_ix(6);
	expected_desc_ix << 0, 2, 4, 1, 3, 5;

	igl::sort(X, 1, true, Y, IX);
	REQUIRE (IX == expected_asc_ix);
	igl::sort(X, 1, false, Y, IX);
	REQUIRE (IX == expected_desc_ix);

	igl::sort_new(X, 1, true, Y, IX);
	REQUIRE (IX == expected_asc_ix);
	igl::sort_new(X, 1, false, Y, IX);
	REQUIRE (IX == expected_desc_ix);

	Eigen::RowVectorXi X_row = X.transpose();
	Eigen::RowVectorXi Y_row, IX_row;
	igl::sort(X_row, 2, true, Y_row, IX_row);
	REQUIRE (IX_row == expected_asc_ix.transpose());
	igl::sort(X_row, 2, false, Y_row, IX_row);
	REQUIRE (IX_row == expected_desc_ix.transpose());

	igl::sort_new(X_row, 2, true, Y_row, IX_row);
	REQUIRE (IX_row == expected_asc_ix.transpose());
	igl::sort_new(X_row, 2, false, Y_row, IX_row);
	REQUIRE (IX_row == expected_desc_ix.transpose());

	Eigen::MatrixXi R(6, 2);
	R << 2, 0,
	     1, 0,
	     2, 0,
	     1, 0,
	     2, 0,
	     1, 0;
	Eigen::MatrixXi SR;
	igl::sortrows(R, true, SR, IX);
	REQUIRE (IX == expected_asc_ix);
	igl::sortrows(R, false, SR, IX);
	REQUIRE (IX == expected_desc_ix);
}
