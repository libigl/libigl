#include <test_common.h>
#include <igl/tet_tet_adjacency.h>
#include <igl/readMESH.h>
#include <iostream>


class tet_tet_adjacency : public ::testing::TestWithParam<std::string> {};

TEST_P(tet_tet_adjacency, dot)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F, T, TT,TTi;
  // Load example mesh: GetParam() will be name of mesh file
  igl::readMESH(test_common::data_path(GetParam()), V, T, F);
  igl::tet_tet_adjacency(T, TT, TTi);
  ASSERT_EQ(T.rows(), TT.rows());
  ASSERT_EQ(T.rows(), TTi.rows());
  ASSERT_EQ(T.cols(),TT.cols());
  ASSERT_EQ(T.cols(),TTi.cols());
  for(int t = 0;t<T.rows();t++)
  {
    for(int c = 0; c<4 ;c++)
    {
      if(TT(t, c) >= 0)
      {
        ASSERT_LT(TT(t, c), T.rows());
        ASSERT_GE(TTi(t, c), 0);
        ASSERT_LT(TTi(t, c), 4);
        ASSERT_EQ(TT(TT(t, c), TTi(t,c)) , t);
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P
(
  tet_meshes,
  tet_tet_adjacency,
  ::testing::ValuesIn(test_common::tet_meshes()),
  test_common::string_test_name
);
