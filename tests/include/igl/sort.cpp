#include <test_common.h>
#include <igl/sort.h>
#include <igl/STR.h>
#include <tuple>

namespace sort
{
  typedef std::tuple<int/*n*/,int/*m*/,int/*dim*/,bool/*ascending*/>
    NMDimAscending;
  inline std::string NMDimAscending_test_name(
    const ::testing::TestParamInfo<NMDimAscending>& info)
  {
    return STR(
      std::get<0>(info.param)<<"x"<<
      std::get<1>(info.param)<<"_"<<
      "dim_"<<std::get<2>(info.param)<<"_"<<
      "ascending_"<<(std::get<3>(info.param)?"true":"false"));
  };
}
class SortTest : public ::testing::TestWithParam<sort::NMDimAscending> {};

TEST_P(SortTest,random)
{
  const int n = std::get<0>(GetParam());
  const int m = std::get<1>(GetParam());
  const int dim = std::get<2>(GetParam());
  const bool ascending = std::get<3>(GetParam());
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(n,m);
  // sort ascending
  Eigen::MatrixXd Y;
  Eigen::MatrixXi IX;
  igl::sort(X,dim,ascending,Y,IX);
  ASSERT_EQ(X.rows(),Y.rows());
  ASSERT_EQ(X.cols(),Y.cols());
  ASSERT_EQ(X.rows(),IX.rows());
  ASSERT_EQ(X.cols(),IX.cols());
  for(int i = 0;i<n;i++)
  {
    for(int j = 0;j<m;j++)
    {
      ASSERT_EQ(Y(i,j),X(dim==1?IX(i,j):i,dim==2?IX(i,j):j));
    }
  }
  for(int i = (dim==1?1:0);i<n;i++) 
  {
    for(int j = (dim==2?1:0);j<m;j++)
    {
      if(ascending)
      {
        ASSERT_LE(Y(i-(dim==1?1:0),j-(dim==2?1:0)),Y(i,j));
      }else
      {
        ASSERT_GE(Y(i-(dim==1?1:0),j-(dim==2?1:0)),Y(i,j));
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P
(
  suite,
  SortTest,
  ::testing::ValuesIn<std::vector<sort::NMDimAscending> >(
    std::vector<sort::NMDimAscending> {
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
    }),
  sort::NMDimAscending_test_name
);
