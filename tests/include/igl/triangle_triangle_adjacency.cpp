
#include <test_common.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>

class triangle_triangle_adjacency : public ::testing::TestWithParam<std::string> {};

TEST_P(triangle_triangle_adjacency, dot)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F,TT,TTi;
  // Load example mesh: GetParam() will be name of mesh file
  test_common::load_mesh(GetParam(), V, F);
  igl::triangle_triangle_adjacency(F,TT,TTi);
  ASSERT_EQ(F.rows(),TT.rows());
  ASSERT_EQ(F.rows(),TTi.rows());
  ASSERT_EQ(F.cols(),TT.cols());
  ASSERT_EQ(F.cols(),TTi.cols());
  for(int f = 0;f<F.rows();f++)
  {
    for(int c = 0;c<3;c++)
    {
      if(TT(f,c) >= 0)
      {
        ASSERT_LT(TT(f,c) , F.rows());
        ASSERT_GE(TTi(f,c) , 0);
        ASSERT_LT(TTi(f,c) , 3);
        ASSERT_EQ( TT(TT(f,c),TTi(f,c)) , f);
      }
    }
  }
  // ASSERT_EQ(a,b);
  // ASSERT_TRUE(a==b);
  // ASSERT_NEAR(a,b,1e-15)
  // ASSERT_LT(a,1e-12);
}

INSTANTIATE_TEST_CASE_P
(
  manifold_meshes,
  triangle_triangle_adjacency,
  ::testing::ValuesIn(test_common::manifold_meshes()),
  test_common::string_test_name
);
