
#include <test_common.h>
#include <igl/is_symmetric.h>

TEST(is_symmetric, sparse)
{
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,4,5,6,7,8,9;
    Eigen::SparseMatrix<double> S = M.sparseView();
    ASSERT_FALSE(igl::is_symmetric(S));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,2,4,5,3,5,6;
    Eigen::SparseMatrix<double> S = M.sparseView();
    ASSERT_FALSE(igl::is_symmetric(S));
  }
}

TEST(is_symmetric, dense)
{
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,4,5,6,7,8,9;
    ASSERT_FALSE(igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,
       2,4,5,
       3,5,6;
    ASSERT_FALSE(igl::is_symmetric(M));
  }
}
