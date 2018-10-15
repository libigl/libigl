#include <test_common.h>
#include <igl/edge_flaps.h>

class edge_flaps : public ::testing::TestWithParam<std::string> {};

TEST_P(edge_flaps, verify)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  test_common::load_mesh(GetParam(), V, F);

  Eigen::MatrixXi efE,efEF,efEI;
  Eigen::VectorXi efEMAP;
  igl::edge_flaps(F,efE,efEMAP,efEF,efEI);
  ASSERT_EQ(efE.rows(),efEF.rows());
  ASSERT_EQ(efE.cols(),2);
  ASSERT_EQ(efE.cols(),efEF.cols());
  // for each edge, make sure edge appears in face
  for(int e = 0;e<efE.rows();e++)
  {
    for(int fe = 0;fe<2;fe++)
    {
      const int f = efEF(e,fe);
      // index of corner
      const int c = efEI(e,fe);
      ASSERT_TRUE(f<F.rows());
      // only check if not on boundary
      if(f >= 0)
      {
        EXPECT_TRUE( 
        // Either efE(e,[1 2]) = [i,j] appears after vertex c of face f
          ((efE(e,0) == F(f,(c+1)%3)) && (efE(e,1) == F(f,(c+2)%3))) ||
        // Or  efE(e,[2 1]) = [j,i] appears after vertex c of face f
          ((efE(e,1) == F(f,(c+1)%3)) && (efE(e,0) == F(f,(c+2)%3))));
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P
(
 all_meshes,
 edge_flaps,
 ::testing::ValuesIn(test_common::all_meshes()),
 test_common::string_test_name
);
